"""Build a joined metadata table for the Shao et al. 2019 infant gut metagenomes.

Joins the ENA run report (technical fields, FTP paths) with the clinical
metadata published as Supplementary Table 1 of the Nature paper.

Inputs
------
ena_link.tsv    ENA filereport, result=read_run, study ERP115334
supp_table.xlsx Supplementary Table 1 (header row is the 3rd row)

Output
------
A TSV with one row per metagenomic run, carrying delivery mode, time point
and FTP download paths.
"""

import argparse
from pathlib import Path

import pandas as pd

SUPP_HEADER_ROW = 2  # 0-indexed; rows 0-1 are the table caption and a blank


def load_ena(path: Path) -> pd.DataFrame:
    """Read the ENA run report and keep metagenomic runs only."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    total = len(df)
    df = df[df["library_source"] == "METAGENOMIC"].copy()
    print(f"ENA report: {total} runs, {len(df)} metagenomic")
    return df


def load_supp(path: Path) -> pd.DataFrame:
    """Read Supplementary Table 1, skipping the caption rows."""
    df = pd.read_excel(path, header=SUPP_HEADER_ROW, dtype=str)
    df.columns = [str(c).strip() for c in df.columns]
    print(f"Supplementary table: {len(df)} rows, {df.shape[1]} columns")
    print(f"Columns: {list(df.columns)}")
    return df


def join(ena: pd.DataFrame, supp: pd.DataFrame) -> pd.DataFrame:
    """Join on the ENA secondary sample accession (ERS...)."""
    merged = ena.merge(
        supp,
        how="left",
        left_on="secondary_sample_accession",
        right_on="Accession",
        validate="one_to_one",
    )
    matched = merged["Accession"].notna().sum()
    print(f"Join: {matched}/{len(merged)} runs matched to clinical metadata")
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ena", type=Path, required=True)
    parser.add_argument("--supp", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    ena = load_ena(args.ena)
    supp = load_supp(args.supp)
    merged = join(ena, supp)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out} ({len(merged)} rows)")


if __name__ == "__main__":
    main()

