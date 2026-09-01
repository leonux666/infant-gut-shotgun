"""Build a long-format species-level abundance table from MetaPhlAn output.

Input
-----
1. MetaPhlAn combined report produced by nf-core/taxprofiler
   (results/metaphlan/*_combined_reports.txt). Line 1 is a database-version
   comment; the header is on line 2. Rows are cumulative across taxonomic
   ranks, so a single rank must be selected before any statistic is computed.
2. Sample metadata (data/selected_samples.tsv), keyed on run_accession.

Output
------
Long-format TSV: one row per sample-by-species pair, carrying Delivery_mode
and Time_point. Long format is chosen because every downstream consumer
(grouped statistics, plotting) requires it; a wide matrix would have to be
melted anyway.

Design notes
------------
- The combined report is used rather than the taxpasta merged table. The
  taxpasta table carries only taxonomy_id, with no lineage string, and
  contains 399 of the 710 rows present here. Three candidate explanations for
  the missing rows were tested and all were ruled out (SGB-level rows alone,
  unnamed clades alone, duplicate-ID collapse), so the selection rule behind
  that table is undocumented and unverified. An input whose filtering rule
  cannot be reproduced is not a sound analysis entry point.
- Eukaryotic clades are retained. MetaPhlAn relative abundances are
  normalised across all kingdoms, so dropping one kingdom would require
  renormalisation, which applies an unequal adjustment to each sample.
  Retained, the species-level matrix sums to 100 per sample as delivered.
- Unnamed clades (e.g. GGB3109_SGB4121) are retained. They carry real
  abundance and occupy a real share of the denominator; they are simply not
  interpretable by name.

Usage
-----
python src/build_species_table.py
python src/build_species_table.py --profile <path> --metadata <path> --output <path>
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import pandas as pd

# Species level is the 7th rank (k|p|c|o|f|g|s), i.e. 6 pipe separators.
SPECIES_PIPE_COUNT = 6

# Line 1 of the combined report is a "#<database version>" comment; the header
# is on line 2. Confirmed empirically, not assumed from documentation.
COMMENT_LINES = 1

# MetaPhlAn relative abundances are percentages summing to 100 per sample
# within a single rank. Tolerance covers float rounding in the source file.
ABUNDANCE_SUM = 100.0
ABUNDANCE_TOL = 0.01

DEFAULT_PROFILE_GLOB = "results/metaphlan/*_combined_reports.txt"
DEFAULT_METADATA = "data/selected_samples.tsv"
DEFAULT_OUTPUT = "results/analysis/species_abundance.tsv"

METADATA_COLS = ["run_accession", "Delivery_mode", "Time_point"]


def resolve_profile(pattern: str) -> Path:
    """Resolve the combined report by glob.

    The filename embeds db_name (e.g. mpa_vJun23_202403), which changes with
    the database version, so it must not be hard-coded.
    """
    matches = sorted(glob.glob(pattern))
    if len(matches) != 1:
        sys.exit(
            f"ERROR: expected exactly 1 file matching {pattern!r}, found {len(matches)}"
        )
    path = Path(matches[0])
    print(f"Profile resolved: {path}")
    return path


def load_profile(path: Path) -> pd.DataFrame:
    """Read the combined report, skipping the database-version comment line."""
    df = pd.read_csv(path, sep="\t", skiprows=COMMENT_LINES, dtype=str)
    print(f"Rows read from profile: {len(df)}")
    return df


def select_species(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only species-level rows.

    Rows in a MetaPhlAn profile are cumulative across ranks: k__Bacteria
    equals the sum of everything beneath it. Selecting one rank first is
    mandatory, otherwise abundances are counted many times over.
    """
    clade_col = df.columns[0]
    is_species = df[clade_col].str.count(r"\|") == SPECIES_PIPE_COUNT
    out = df.loc[is_species].copy()
    print(f"Species-level rows retained: {len(out)} of {len(df)}")
    if out.empty:
        sys.exit("ERROR: no species-level rows found; check the lineage format")
    return out


def parse_sample_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Rename abundance columns to bare run accessions.

    Column names arrive as
    '<sample>_<run>_<db_name>.metaphlan', and sample and run were set to the
    same value in the samplesheet, so the accession appears twice. Taking the
    leading underscore-delimited field recovers the accession.
    """
    clade_col = df.columns[0]
    renamed = {clade_col: "clade_name"}
    for col in df.columns[1:]:
        renamed[col] = col.split("_")[0]
    out = df.rename(columns=renamed)

    sample_cols = list(out.columns[1:])
    print(f"Sample columns parsed: {len(sample_cols)}")
    if len(set(sample_cols)) != len(sample_cols):
        sys.exit("ERROR: sample column names are not unique after parsing")
    return out


def to_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """Cast abundances to float and verify each sample sums to 100.

    This is the strongest self-check in the script. A rank-selection error
    inflates the sum (multiple ranks stacked); a parsing error deflates it
    (rows lost). One number catches both failure modes.
    """
    out = df.copy()
    sample_cols = list(out.columns[1:])
    out[sample_cols] = out[sample_cols].astype(float)

    sums = out[sample_cols].sum()
    off = sums[(sums - ABUNDANCE_SUM).abs() > ABUNDANCE_TOL]
    print(f"Samples summing to {ABUNDANCE_SUM} +/- {ABUNDANCE_TOL}: "
          f"{len(sample_cols) - len(off)} of {len(sample_cols)}")
    if not off.empty:
        for name, value in off.items():
            print(f"  {name}: {value}")
        sys.exit("ERROR: abundance sums are off; rank selection or parsing is wrong")
    return out


def melt_and_join(df: pd.DataFrame, metadata_path: Path) -> pd.DataFrame:
    """Reshape to long format and attach Delivery_mode and Time_point."""
    long_df = df.melt(
        id_vars="clade_name",
        var_name="run_accession",
        value_name="relative_abundance",
    )
    print(f"Long-format rows: {len(long_df)}")

    # Species name is the terminal lineage field, with its 's__' prefix removed.
    long_df["species"] = (
        long_df["clade_name"].str.split("|").str[-1].str.replace("s__", "", regex=False)
    )

    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)
    missing = [c for c in METADATA_COLS if c not in meta.columns]
    if missing:
        sys.exit(f"ERROR: metadata is missing columns: {missing}")
    meta = meta[METADATA_COLS]
    print(f"Metadata rows read: {len(meta)}")

    merged = long_df.merge(meta, on="run_accession", how="left", validate="many_to_one")

    unmatched = merged["Delivery_mode"].isna().sum()
    print(f"Rows with no metadata match: {unmatched}")
    if unmatched:
        sys.exit("ERROR: some samples have no metadata; check run_accession values")

    return merged[
        ["run_accession", "Delivery_mode", "Time_point", "species", "relative_abundance"]
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--profile", default=DEFAULT_PROFILE_GLOB,
                        help=f"glob for the combined report (default: {DEFAULT_PROFILE_GLOB})")
    parser.add_argument("--metadata", default=DEFAULT_METADATA,
                        help=f"sample metadata TSV (default: {DEFAULT_METADATA})")
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help=f"output TSV (default: {DEFAULT_OUTPUT})")
    args = parser.parse_args()

    profile_path = resolve_profile(args.profile)
    metadata_path = Path(args.metadata)
    if not metadata_path.exists():
        sys.exit(f"ERROR: metadata not found: {metadata_path}")

    df = load_profile(profile_path)
    df = select_species(df)
    df = parse_sample_ids(df)
    df = to_numeric(df)
    result = melt_and_join(df, metadata_path)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {len(result)} rows to {output_path}")


if __name__ == "__main__":
    main()
