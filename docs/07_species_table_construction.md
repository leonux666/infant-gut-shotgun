# 07 — Species table construction

Builds the species-level abundance matrix that downstream comparison and
plotting consume. This is the first step in the project that leaves the
pipeline and enters analysis: everything before it produced taxonomic
profiles, and everything after it interprets them.

Script: `src/build_species_table.py`
Output: `results/analysis/species_abundance.tsv` (not in git; rebuildable)

## Input file structure

Two facts about the MetaPhlAn combined report have to be handled explicitly,
and neither is visible from the filename or the nf-core documentation.

**The header is on line 2.** Line 1 is a comment naming the database version
(`#mpa_vJun23_CHOCOPhlAnSGB_202403`). Reading with a default `header=0`
silently produces a one-column frame, because the comment line contains no
tab. The per-sample profiles under `metaphlan/<db_name>/` have a different
comment block again — five lines, the fifth being the header — so the two
file types cannot share a reader.

**Rows are cumulative across ranks.** The row for `k__Bacteria` equals the sum
of every bacterial row beneath it. A statistic computed over the unfiltered
table counts the same abundance once per rank. One rank must be selected
before anything else happens.

Rank depth is encoded in the pipe count of the lineage string. Observed
distribution across the 710 data rows:

| Depth | Rank | Rows |
|---|---|---|
| 1 | kingdom | 2 |
| 2 | phylum | 6 |
| 3 | class | 18 |
| 4 | order | 29 |
| 5 | family | 49 |
| 6 | genus | 112 |
| 7 | species | 239 |
| 8 | SGB | 255 |

Abundance column names arrive as `<sample>_<run>_<db_name>.metaphlan`. Sample
and run were set to the same value in the samplesheet, so the run accession
appears twice; the leading underscore-delimited field recovers it. The
`db_name` is embedded in both the column names and the filename, so neither
may be hard-coded — the script resolves the file by glob.

## Choosing the analysis entry point

The pipeline produces two candidate merged tables. The choice between them
was not obvious and is recorded here because it changes what downstream code
can see.

| | `taxpasta/*.tsv` | `metaphlan/*_combined_reports.txt` |
|---|---|---|
| First column | `taxonomy_id` only | `clade_name`, full lineage |
| Data rows | 399 | 710 |
| Rank recoverable | no | yes, from pipe count |

The taxpasta table carries no lineage string. Rank cannot be recovered from a
bare NCBI taxonomy ID without an external mapping, and rank selection is
mandatory here, so that table cannot serve as the entry point unaided.

The row count gap was investigated before drawing any conclusion, because a
table holding 399 of 710 rows is applying a filter, and an unexplained filter
on the input is worse than a missing column. Three candidate explanations
were tested; all three failed.

**Hypothesis 1 — the table excludes SGB-level rows.**
Symptom: 710 minus the 255 depth-8 rows leaves 455, not 399.
Verdict: rejected, off by 56.

**Hypothesis 2 — the table additionally excludes unnamed clades.**
MetaPhlAn assigns placeholder names (`c__CFGB1292`, `g__GGB3109`,
`s__GGB3109_SGB4121`) to clades with no formal taxonomy. Counting these
under two different pattern definitions gave 51 and 36, neither of which
closes a gap of 56.
Verdict: rejected.

**Hypothesis 3 — rows sharing a taxonomy ID were collapsed.**
`s__Candida_albicans` and `t__EUK5476` share NCBI taxid 5476, so a table keyed
on taxonomy ID could merge such pairs. Checking the taxpasta ID column for
duplicates returned zero, meaning no collapse occurred — but also that a
one-to-one keying would not have lost rows this way.
Verdict: rejected.

A fourth check tested whether unnamed clades simply lack an NCBI ID and are
dropped for that reason. In one per-sample profile, 74 non-SGB rows were
present and zero had an empty terminal taxonomy ID, so MetaPhlAn assigns IDs
to placeholder clades as well.

**Conclusion.** The selection rule behind the taxpasta table could not be
reproduced from the data. That is itself the finding: an input whose
filtering cannot be explained is not a sound analysis entry point, regardless
of whether the filtering is correct. The combined report is used instead —
it carries full lineage, its row count agrees with the per-sample profiles,
and it applies no transformation that this repository cannot account for.
The taxpasta table is retained in `results/` for cross-checking and is not
read by any script here.

## Scope decisions

**Species rank only.** Genus and phylum tables are equally easy to produce
from the same input by changing one constant, but each additional rank
multiplies the downstream comparisons without adding a distinct question.
Species is the rank at which the birth-mode literature reports its findings,
so it is the rank that can be checked against published work.

**Eukaryotic clades retained.** The profile contains one eukaryotic lineage,
`Candida albicans`, present at species rank. MetaPhlAn normalises relative
abundance across all kingdoms, so removing a kingdom requires renormalising
to 100, and that renormalisation applies a different adjustment to every
sample depending on how much was removed there. Retained, each sample's
species vector sums to 100 exactly as delivered and no transformation is
introduced. The cost is that the matrix describes the microbiome rather than
the bacterial community, and downstream wording must say so.

**Unnamed clades retained.** Roughly a dozen species-rank rows carry
placeholder names of the form `GGB3109_SGB4121`. They hold real abundance and
occupy a real share of the denominator, so dropping them would again require
renormalisation. They are kept in the matrix and must be labelled as
unnamed wherever taxa are listed or plotted — they are not interpretable by
name and should not be presented as though they were.

## Verification

The script prints one checkable count per stage. These are not decoration;
each has a specific failure it is positioned to catch.

| Count | Expected | Catches |
|---|---|---|
| Rows read from profile | 710 | wrong `skiprows` |
| Species-level rows retained | 239 | wrong depth predicate |
| Sample columns parsed | 24 | column-name parsing error |
| Samples summing to 100 | 24 of 24 | see below |
| Long-format rows | 5736 | rows lost in reshape |
| Rows with no metadata match | 0 | accession mismatch |

The per-sample sum is the strongest of these because it catches two opposite
failures with one number. Stacking more than one rank inflates the sum above
100; losing rows during parsing deflates it below. Any deviation beyond 0.01
exits non-zero rather than writing a plausible-looking but wrong table.

## Reproduce

    conda activate nfcore
    python src/build_species_table.py

Defaults resolve the combined report by glob under `results/metaphlan/`, read
`data/selected_samples.tsv`, and write `results/analysis/species_abundance.tsv`.
All three are overridable with `--profile`, `--metadata` and `--output`.

Output is long-format TSV with five columns: `run_accession`,
`Delivery_mode`, `Time_point`, `species`, `relative_abundance`. Long format
is used because grouped statistics and plotting both require it; a wide
matrix would have to be reshaped anyway.

## Not yet done

- No statistical comparison. With three samples per delivery-mode by
  time-point cell, the design supports description and ordination but not a
  differential-abundance test that would survive scrutiny. Any comparison
  made from this table must be presented as descriptive.
- Genus and phylum tables are not produced. If they are wanted later, the
  change is the `SPECIES_PIPE_COUNT` constant plus the expected-count checks
  that depend on it.
- Prevalence and abundance filtering is not applied. The table holds every
  species MetaPhlAn reported, including those detected in a single sample.
- No plots. Downstream visualisation lives in `notebooks/` and is not yet
  written.
