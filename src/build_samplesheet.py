"""Write the nf-core/taxprofiler input samplesheet.

Input : data/selected_samples.tsv  (24 rows, from select_samples.py)
Output: data/samplesheet.csv       (24 rows, six columns)

Column layout is fixed by the pipeline (taxprofiler 2.0.1):

    sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta

'sample' identifies the biological sample and 'run_accession' the sequencing
run. The pipeline concatenates runs sharing a sample identifier when
--perform_runmerging is set. Each of these samples was sequenced once, so the
two columns carry the same value and run merging stays off.

'fasta' is for pre-assembled long-read input and is left empty.

Paths are absolute and point at the real files on scratch rather than at the
data/fastq symlink, since Nextflow stages inputs from wherever the path
resolves and an unresolved link is one more thing that can go wrong.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Fixed by the pipeline's input schema; order matters.
COLUMNS = [
    "sample",
    "run_accession",
    "instrument_platform",
    "fastq_1",
    "fastq_2",
    "fasta",
]

# ENA reports the two mates of a paired-end run as one semicolon-separated
# field, R1 first.
FIELD_SEP = ";"

# All 1,679 runs in ERP115334 are Illumina; taxprofiler expects this spelling.
PLATFORM = "ILLUMINA"


def load_runs(path: Path) -> pd.DataFrame:
    """Read the selected-sample table."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    print(f"[load] {len(df)} runs")
    return df


def build_rows(df: pd.DataFrame, fastq_dir: Path) -> pd.DataFrame:
    """Map each run to one samplesheet row with absolute FASTQ paths."""
    rows = []
    for _, run in df.iterrows():
        urls = run["fastq_ftp"].split(FIELD_SEP)
        if len(urls) != 2:
            raise ValueError(
                f"{run['run_accession']}: expected 2 FASTQ files, got {len(urls)}"
            )
        r1, r2 = (fastq_dir / Path(u).name for u in urls)
        for mate in (r1, r2):
            if not mate.exists():
                raise FileNotFoundError(mate)
        rows.append(
            {
                "sample": run["run_accession"],
                "run_accession": run["run_accession"],
                "instrument_platform": PLATFORM,
                "fastq_1": str(r1),
                "fastq_2": str(r2),
                "fasta": "",
            }
        )

    sheet = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[build] {len(sheet)} rows, {2 * len(sheet)} FASTQ files verified")
    return sheet


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--fastq-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    df = load_runs(args.input)
    sheet = build_rows(df, args.fastq_dir.resolve())
    sheet.to_csv(args.output, index=False)
    print(f"[write] {args.output} ({len(sheet)} rows)")


if __name__ == "__main__":
    sys.exit(main())
