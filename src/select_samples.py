"""Select a balanced subset of infant gut metagenome samples for profiling.

Input : data/metadata_joined.tsv  (1679 rows, ENA run report joined with
        Shao et al. 2019 Supplementary Table 1)
Output: data/selected_samples.tsv (24 rows, one row per infant)

Design
------
Cross-sectional, not longitudinal: every selected infant appears exactly once.
Cells are fully crossed and equally sized so that delivery mode cannot be
confounded with sampling time point.

    2 delivery modes x 4 time points x 3 samples = 24

Infant in-hospital antibiotic exposure is held constant rather than balanced:
only unexposed infants are eligible. Exposure runs at 5-11% across the cells,
so balancing 1-of-3 would have oversampled it several fold, and neonatal
antibiotics perturb colonisation in the same direction as caesarean birth.
At three samples per cell the covariate cannot be estimated anyway, so
restriction buys lower within-cell variance at almost no cost in
representativeness.

Why 24: 1.22 GB mean paired-end size per run, so the subset is roughly 30 GB.
This is a compute and storage budget, not a power calculation.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Supplementary Table 1 records four scheduled visits plus maternal samples and
# a handful of off-schedule infant days. Only the scheduled infant visits are
# used so that the design stays balanced.
TIME_POINTS = ["4", "7", "21", "Infancy"]
DELIVERY_MODES = ["Caesarean", "Vaginal"]
N_PER_CELL = 3

# Abx_Baby_in_hospital is 94% complete (1386/1469) but mixes 'No' and 'no'.
# Abx_mother_labour_IAP was rejected instead: 16% complete, and zero 'Yes'
# in the caesarean arm, so it is not comparable across delivery modes.
ABX_COL = "Abx_Baby_in_hospital"
ABX_KEEP = "no"

# Stratified sampling draws each cell independently, so one infant can be
# picked twice via two different time points. With ~500 infants contributing
# ~1400 records this happens often, not rarely: seeds 42-52 all collided
# under the earlier design. Reseed and redraw rather than swapping a record,
# which would break the stratification.
MAX_SEED_ATTEMPTS = 100


def load_metadata(path: Path) -> pd.DataFrame:
    """Read the joined metadata table, keeping every field as a string."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    print(f"[load] {len(df)} rows, {df.shape[1]} columns")
    return df


def filter_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """Restrict to scheduled infant visits without antibiotic exposure."""
    df = df[df["Time_point"].isin(TIME_POINTS)].copy()
    print(f"[filter] {len(df)} rows at scheduled infant time points")

    df["abx"] = df[ABX_COL].str.strip().str.lower()
    df = df[df["abx"] == ABX_KEEP]
    print(f"[filter] {len(df)} rows without in-hospital antibiotics")

    return df


def draw_subset(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Draw one balanced subset; may contain repeated individuals."""
    picks = []
    for mode in DELIVERY_MODES:
        for tp in TIME_POINTS:
            pool = df[
                (df["Delivery_mode"] == mode) & (df["Time_point"] == tp)
            ]
            if len(pool) < N_PER_CELL:
                raise ValueError(
                    f"cell {mode}/{tp} holds {len(pool)} rows, "
                    f"need {N_PER_CELL}"
                )
            picks.append(pool.sample(n=N_PER_CELL, random_state=seed))
    return pd.concat(picks)


def select(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Redraw until every selected infant is distinct."""
    for attempt in range(MAX_SEED_ATTEMPTS):
        subset = draw_subset(df, seed + attempt)
        n_individuals = subset["Individual"].nunique()
        if n_individuals == len(subset):
            print(f"[select] {len(subset)} rows, {n_individuals} individuals, "
                  f"seed {seed + attempt}")
            return subset
        print(f"[select] seed {seed + attempt} repeated an individual "
              f"({n_individuals} of {len(subset)}), redrawing")

    raise RuntimeError(
        f"no distinct-individual draw within {MAX_SEED_ATTEMPTS} seeds"
    )


def report(subset: pd.DataFrame) -> None:
    """Print the realised design so the balance can be checked by eye."""
    counts = subset.groupby(["Delivery_mode", "Time_point"]).size()
    print("\n[design]")
    print(counts.to_string())

    # ENA reports paired-end runs as 'R1_bytes;R2_bytes' in a single field.
    total_bytes = (
        subset["fastq_bytes"]
        .str.split(";")
        .explode()
        .astype(float)
        .sum()
    )
    print(f"\n[design] total download: {total_bytes / 1e9:.1f} GB")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    df = load_metadata(args.input)
    df = filter_candidates(df)
    subset = select(df, args.seed)
    report(subset)

    subset = subset.drop(columns=["abx"])
    subset.to_csv(args.output, sep="\t", index=False)
    print(f"\n[write] {args.output} ({len(subset)} rows)")


if __name__ == "__main__":
    sys.exit(main())
