"""Download the selected FASTQ files from the ENA FTP mirror.

Input : data/selected_samples.tsv  (24 rows, from select_samples.py)
Output: 48 gzipped FASTQ files in the target directory

Files are fetched one at a time with wget. ENA serves a public archive, so
there is no reason to open parallel connections for a set this small.

Verification is by byte count against the fastq_bytes field, which catches
truncated transfers. Silent corruption at matching size is rare over HTTP and
would surface later as a read error during quality trimming.

Already-complete files are skipped, so the script can be rerun after an
interrupted job without refetching what is already on disk.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

# ENA reports the two mates of a paired-end run as one semicolon-separated
# field, in both fastq_ftp and fastq_bytes, in the same order.
FIELD_SEP = ";"

# fastq_ftp omits the scheme. HTTPS is used rather than FTP because outbound
# FTP is filtered on many clusters.
URL_SCHEME = "https://"

# wget flags: -c resumes a partial file, -q suppresses the progress bar that
# would otherwise fill the SLURM log, --tries retries a dropped connection.
WGET_FLAGS = ["-c", "-q", "--tries=3"]


def load_targets(path: Path) -> list[tuple[str, int]]:
    """Return one (url, expected_bytes) pair per FASTQ file."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    print(f"[load] {len(df)} runs")

    targets: list[tuple[str, int]] = []
    for _, row in df.iterrows():
        urls = row["fastq_ftp"].split(FIELD_SEP)
        sizes = row["fastq_bytes"].split(FIELD_SEP)
        if len(urls) != len(sizes):
            raise ValueError(
                f"{row['run_accession']}: {len(urls)} urls, {len(sizes)} sizes"
            )
        targets.extend(zip(urls, (int(s) for s in sizes)))

    total_gb = sum(size for _, size in targets) / 1e9
    print(f"[load] {len(targets)} files, {total_gb:.1f} GB")
    return targets


def fetch(url: str, expected: int, outdir: Path) -> str:
    """Download one file unless it is already present at full size."""
    dest = outdir / Path(url).name

    if dest.exists() and dest.stat().st_size == expected:
        return "skipped"

    subprocess.run(
        ["wget", *WGET_FLAGS, "-O", str(dest), URL_SCHEME + url],
        check=True,
    )

    actual = dest.stat().st_size
    if actual != expected:
        raise ValueError(
            f"{dest.name}: got {actual} bytes, expected {expected}"
        )
    return "fetched"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    targets = load_targets(args.input)

    counts = {"fetched": 0, "skipped": 0}
    for i, (url, expected) in enumerate(targets, start=1):
        status = fetch(url, expected, args.outdir)
        counts[status] += 1
        print(f"[{i}/{len(targets)}] {Path(url).name} {status}", flush=True)

    print(f"\n[done] {counts['fetched']} fetched, {counts['skipped']} skipped")


if __name__ == "__main__":
    sys.exit(main())
