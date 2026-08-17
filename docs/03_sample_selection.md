# 03 — Sample selection

Selects 24 runs from the 1,679 infant gut metagenomes in ERP115334 for
taxonomic and functional profiling.

- Script: `src/select_samples.py`
- Input: `data/metadata_joined.tsv` (see `02_metadata_assembly.md`)
- Output: `data/selected_samples.tsv`

```bash
python src/select_samples.py \
  --input data/metadata_joined.tsv \
  --output data/selected_samples.tsv
```

The seed defaults to 42 and is not passed on the command line, so the same
command reproduces the same 24 runs.

---

## Why 24

Paired-end runs average 1.22 GB, so 24 runs is roughly 30 GB — small enough
for scratch storage and for a full profiling run to finish within the project
timeline. This is a compute budget, not a power calculation, and the analysis
is framed accordingly.

---

## Design

Cross-sectional. Delivery mode and time point are fully crossed with equal
cell sizes:

|           | Day 4 | Day 7 | Day 21 | Infancy |
|-----------|-------|-------|--------|---------|
| Caesarean | 3     | 3     | 3      | 3       |
| Vaginal   | 3     | 3     | 3      | 3       |

Equal cells make it structurally impossible for delivery mode to be
confounded with sampling time.

### One sample per infant

The cohort is longitudinal: the same infant appears at several time points, so
the 1,679 runs come from roughly 500 infants. Drawing each cell independently
can therefore pick one infant twice.

A half-overlapping set is the worst outcome — neither a paired design nor
independent observations, and awkward to analyse either way. The script
asserts that the number of distinct `Individual` values equals the number of
selected rows, and reseeds and redraws when it does not.

Collisions are common rather than rare, which is why the assertion is
necessary rather than defensive: under an earlier variant that drew twice per
cell, seeds 42 through 52 all collided. The current design collides on the
first seed and succeeds on the second.

A longitudinal alternative — 6 infants followed across all four visits — was
considered and rejected: 6 individuals cannot support any inference about
group differences.

---

## Antibiotic exposure: held constant, not balanced

Only infants recorded as unexposed to in-hospital antibiotics are eligible
(`Abx_Baby_in_hospital == "no"`, case-normalised). This removes 130 of the
1,386 resolved records, leaving 1,256 candidates.

Three reasons:

1. **Prevalence.** Exposure runs at 5–11% across the eight cells. Balancing
   one exposed sample into every cell of three would have set it at 33%,
   oversampling a minority condition several fold.
2. **Direction.** Neonatal antibiotics shift colonisation the same way
   caesarean birth does — suppressing *Bifidobacterium* and *Bacteroides*,
   favouring *Enterococcus* and *Klebsiella*. Balanced allocation keeps the
   group means unbiased but blunts the contrast, which is scarce at n=24.
3. **Estimability.** With three samples per cell the covariate cannot be
   estimated, so balancing pays a cost in variance without buying anything.

Restriction narrows the population the results speak to — unexposed infants —
which is the large majority of this cohort.

### Why not maternal intrapartum prophylaxis

`Abx_mother_labour_IAP` was the first candidate covariate and was rejected on
inspection:

| Column                       | Complete (of 1,469) |
|------------------------------|---------------------|
| `Abx_mother_prior_birth`     | 1,469               |
| `Abx_mother_labour_IAP`      | 242                 |
| `Abx_mother_after_hospital`  | 1,468               |
| `Abx_Baby_in_hospital`       | 1,386               |
| `Abx_Baby_after_hospital`    | 316                 |

Beyond the 16% completeness, the caesarean arm contains no `Yes` values at
all, which contradicts routine prophylaxis at caesarean delivery. The field
appears to be recorded only for vaginal births, making it incomparable across
the very groups being contrasted.

`Abx_Baby_in_hospital` is also preferable on substance: it is exposure acting
directly on the infant gut, not exposure inferred to reach the infant through
the mother.

---

## Host read removal is switched off

Supplementary Table 1 reports `WGS_reads_raw` and `WGS_reads_trimmed` per run.
Their ratio bounds how much quality trimming and host depletion together
removed, without running anything.

Across all 1,679 runs: mean retention 86.3%, range 67.1–96.9%, single-peaked,
with 92% of runs between 80% and 95%. There is no second mode, so no batch of
runs lost a large fraction.

Retention does not vary with sampling time:

| Time point | n   | Retention |
|------------|-----|-----------|
| Day 4      | 310 | 0.863     |
| Day 7      | 532 | 0.862     |
| Day 21     | 325 | 0.847     |
| Infancy    | 302 | 0.877     |
| Mother     | 175 | 0.862     |

Host DNA in infant stool declines over the first weeks of life, so a
meaningful host-depletion component should leave a gradient. Instead the range
is three percentage points, it is not monotonic, and adult maternal samples
match neonatal ones exactly.

The losses are therefore dominated by quality trimming, and host reads are too
few to detect by this measure. Host removal is left off in the profiling
pipeline: it requires downloading and indexing a human reference for a step
with almost nothing to remove. Quality trimming is left on — it costs minutes
and produces a QC report.

This evidence is indirect. It bounds the combined loss and shows the loss is
uniform; it cannot separate trimmed bases from host reads.

---

## Excluded records

| Excluded                                | Reason                                   |
|-----------------------------------------|------------------------------------------|
| 708 runs with `library_source=GENOMIC`  | Bacterial isolates, not metagenomes       |
| 175 `Mother` samples                    | Not infants                               |
| 35 off-schedule infant days (6, 8–18)   | Too few per day to balance                |
| 83 records with no antibiotic status    | Eligibility cannot be determined          |
| 130 antibiotic-exposed infants          | Excluded by restriction (above)           |
