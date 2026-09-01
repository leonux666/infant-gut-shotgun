"""Plot species composition and key-genus abundance from the species table.

Input
-----
results/analysis/species_abundance.tsv, the long-format table written by
src/build_species_table.py (run_accession, Delivery_mode, Time_point,
species, relative_abundance).

Output
------
Two PNGs under results/analysis/figures/:

  composition_stacked.png  Per-sample stacked bars, top N species plus Other,
                           grouped by delivery mode and ordered by time point.
  key_genera.png           Per-sample abundance of three genera chosen to
                           carry information in both directions: Bacteroides
                           and Bifidobacterium are expected to be reduced
                           after caesarean birth, Escherichia to be enriched.

Design notes
------------
- Top species are ranked by mean abundance across all samples, not by
  maximum. One sample here is 89% a single species; ranking by maximum would
  let such a sample displace taxa that sit at moderate abundance in many
  samples, which is the more informative pattern.
- Genus abundance is derived by prefix match on the species name, because
  the input table holds species rank only. MetaPhlAn species names are
  formatted Genus_species, so the prefix is reliable, but the matched names
  are printed so the match can be checked rather than trusted.
- Time_point values are strings ('4', '7', '21', 'Infancy') and sort into
  the wrong order lexically ('21' before '4'). The order is pinned explicitly.
- Unnamed clades keep their placeholder names in legends. They are real
  abundance and are not interpretable by name; renaming them would hide that.
- Linear y axis, not log. Abundances span zero to 89% and true zeros are
  common, which a log axis cannot show.
- Jitter in the genus panels is random, not index-derived. An index-derived
  offset places tied values on a deterministic line, which hides how many
  points are stacked at zero; with most samples at zero here, that is the
  part of the figure that most needs to be legible.

Usage
-----
python src/plot_composition.py
python src/plot_composition.py --input <path> --outdir <path> --top-n 10
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # No display on a login node; render straight to file.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_INPUT = "results/analysis/species_abundance.tsv"
DEFAULT_OUTDIR = "results/analysis/figures"
DEFAULT_TOP_N = 10

# Sampling ages are stored as strings and would sort as '21' < '4' < '7'.
TIME_ORDER = ["4", "7", "21", "Infancy"]
DELIVERY_ORDER = ["Vaginal", "Caesarean"]

# Two genera expected to be depleted after caesarean birth and one expected
# to be enriched, so the figure can disagree with the expectation.
KEY_GENERA = ["Bacteroides", "Bifidobacterium", "Escherichia"]

OTHER_LABEL = "Other"
EXPECTED_SUM = 100.0
SUM_TOL = 0.01

# Fixed so the jitter is identical on every run; the figure is a committed
# artefact and should not change when nothing upstream has changed.
JITTER_SEED = 42
JITTER_WIDTH = 0.09


def load_table(path: Path) -> pd.DataFrame:
    """Read the long-format species table and check it is internally whole."""
    df = pd.read_csv(path, sep="\t", dtype={"Time_point": str})
    print(f"Rows read: {len(df)}")

    n_samples = df["run_accession"].nunique()
    n_species = df["species"].nunique()
    print(f"Samples: {n_samples}   Species: {n_species}")

    if len(df) != n_samples * n_species:
        sys.exit(
            f"ERROR: table is not rectangular: {len(df)} rows for "
            f"{n_samples} samples x {n_species} species"
        )

    sums = df.groupby("run_accession")["relative_abundance"].sum()
    off = sums[(sums - EXPECTED_SUM).abs() > SUM_TOL]
    print(f"Samples summing to {EXPECTED_SUM}: {n_samples - len(off)} of {n_samples}")
    if not off.empty:
        sys.exit("ERROR: abundances do not sum to 100; rebuild the species table")

    unexpected = set(df["Time_point"]) - set(TIME_ORDER)
    if unexpected:
        sys.exit(f"ERROR: unexpected Time_point values: {sorted(unexpected)}")

    return df


def order_samples(df: pd.DataFrame) -> pd.DataFrame:
    """Order samples by delivery mode, then time point, then accession.

    Both variables are readable off the x axis this way: delivery mode as two
    blocks, time point as a gradient within each block. Returns the ordered
    key frame so callers can reuse the group sizes rather than recompute them.
    """
    key = df[["run_accession", "Delivery_mode", "Time_point"]].drop_duplicates()
    key["_d"] = key["Delivery_mode"].map({v: i for i, v in enumerate(DELIVERY_ORDER)})
    key["_t"] = key["Time_point"].map({v: i for i, v in enumerate(TIME_ORDER)})
    key = key.sort_values(["_d", "_t", "run_accession"]).reset_index(drop=True)

    counts = key["Delivery_mode"].value_counts()
    print("\nSamples per delivery mode:")
    for mode in DELIVERY_ORDER:
        print(f"  {mode:10s} {counts.get(mode, 0)}")

    return key


def top_species(df: pd.DataFrame, top_n: int) -> list[str]:
    """Rank species by mean abundance across all samples."""
    means = df.groupby("species")["relative_abundance"].mean().sort_values(ascending=False)
    selected = means.head(top_n).index.tolist()

    print(f"\nTop {top_n} species by mean abundance:")
    for name in selected:
        print(f"  {means[name]:6.2f}%  {name}")
    print(f"  {means.drop(selected).sum():6.2f}%  ({len(means) - top_n} others)")

    return selected


def plot_stacked(df: pd.DataFrame, top: list[str], key: pd.DataFrame, out: Path) -> None:
    """Stacked bars, one per sample, top species plus a pooled Other."""
    order = key["run_accession"].tolist()
    wide = df.pivot(index="run_accession", columns="species",
                    values="relative_abundance").loc[order]

    plot_df = wide[top].copy()
    plot_df[OTHER_LABEL] = wide.drop(columns=top).sum(axis=1)

    fig, ax = plt.subplots(figsize=(14, 7))
    colours = plt.cm.tab20(range(len(top))).tolist() + [(0.75, 0.75, 0.75, 1.0)]

    bottom = pd.Series(0.0, index=plot_df.index)
    for colour, name in zip(colours, plot_df.columns):
        ax.bar(range(len(plot_df)), plot_df[name], bottom=bottom,
               color=colour, label=name, width=0.85)
        bottom += plot_df[name]

    # Block boundaries come from the ordered key, so an unequal split is
    # placed correctly rather than assumed to be halfway.
    starts: dict[str, int] = {}
    sizes: dict[str, int] = {}
    cursor = 0
    for mode in DELIVERY_ORDER:
        n = int((key["Delivery_mode"] == mode).sum())
        starts[mode], sizes[mode] = cursor, n
        cursor += n

    boundary = starts[DELIVERY_ORDER[1]] - 0.5
    ax.axvline(boundary, color="black", linewidth=1.2)

    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(key["Time_point"].tolist(), fontsize=8)
    ax.set_xlim(-0.6, len(plot_df) - 0.4)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Sampling age (days, or Infancy)")
    ax.set_ylabel("Relative abundance (%)")

    # Group labels sit inside the axes; the title is left to sit above them.
    for mode in DELIVERY_ORDER:
        centre = starts[mode] + sizes[mode] / 2 - 0.5
        ax.text(centre, 97.5, mode, ha="center", va="top", fontsize=12,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="none", alpha=0.85))

    ax.set_title("Species composition by sample", pad=12)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8,
              title="Species", title_fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def plot_key_genera(df: pd.DataFrame, out: Path) -> None:
    """One panel per genus, points per sample, split by delivery mode."""
    rng = np.random.default_rng(JITTER_SEED)
    fig, axes = plt.subplots(1, len(KEY_GENERA), figsize=(14, 4.5))

    print("\nGenus prefix matches:")
    for ax, genus in zip(axes, KEY_GENERA):
        matched = sorted(
            s for s in df["species"].unique() if s.startswith(f"{genus}_")
        )
        print(f"  {genus}: {len(matched)} species")
        for name in matched:
            print(f"      {name}")

        sub = df[df["species"].isin(matched)]
        totals = (
            sub.groupby(["run_accession", "Delivery_mode"])
            ["relative_abundance"].sum().reset_index()
        )

        for i, mode in enumerate(DELIVERY_ORDER):
            vals = totals.loc[totals["Delivery_mode"] == mode, "relative_abundance"]
            x = i + rng.uniform(-JITTER_WIDTH, JITTER_WIDTH, size=len(vals))
            ax.scatter(x, vals, s=45, alpha=0.8,
                       edgecolor="black", linewidth=0.5)
            n_zero = int((vals < 0.01).sum())
            ax.text(i, -0.13, f"{n_zero}/{len(vals)} at 0", ha="center",
                    va="top", fontsize=8, color="dimgray",
                    transform=ax.get_xaxis_transform())

        ax.set_xticks(range(len(DELIVERY_ORDER)))
        ax.set_xticklabels(DELIVERY_ORDER)
        ax.set_xlim(-0.5, len(DELIVERY_ORDER) - 0.5)
        ax.set_title(genus, fontstyle="italic")
        ax.set_ylabel("Relative abundance (%)")
        ax.set_ylim(bottom=-2)

    fig.suptitle("Key genera by delivery mode (n = 12 per group, descriptive only)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N)
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"ERROR: input not found: {input_path}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_table(input_path)
    key = order_samples(df)
    top = top_species(df, args.top_n)

    plot_stacked(df, top, key, outdir / "composition_stacked.png")
    plot_key_genera(df, outdir / "key_genera.png")


if __name__ == "__main__":
    main()
