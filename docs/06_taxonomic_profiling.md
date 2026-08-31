# 06 — Taxonomic profiling

Runs nf-core/taxprofiler 2.0.1 over the 24 selected samples to produce
per-sample taxonomic profiles with MetaPhlAn 4.1.1, plus read-level quality
control. Prerequisites are the samplesheet from `docs/04` and the database and
container environment from `docs/05`.

---

## 1. Cluster configuration

`config/amarel.config` describes the cluster; it carries no analysis parameters.
That separation is deliberate: the config answers "what does this machine look
like" and the job script answers "what are we running this time," so the config
can be reused unchanged by any other pipeline on the same cluster.

Every value in it was measured rather than assumed:

| Setting | Value | Basis |
|---|---|---|
| executor | slurm | Without this Nextflow runs every task locally, on one node |
| queue | main | The default partition, 509 compute nodes, from `sinfo -s` |
| resourceLimits cpus | 24 | Smallest node in main has 32 cores; 24 leaves headroom |
| resourceLimits memory | 120 GB | Smallest node has 192000 MB; 120 leaves headroom |
| resourceLimits time | 24 h | Partition allows 3 days, but a task running a full day here indicates a fault and should fail rather than occupy the queue |
| maxRetries | 2 | — |

Two points that are easy to get wrong:

- taxprofiler 2.0.1 uses the `resourceLimits` list inside a `process` block. The
  older `params.max_cpus` / `max_memory` / `max_time` convention does not exist
  in this version, and setting it has no effect. Confirmed by grepping the
  pipeline's own config files rather than reading the documentation.
- These are ceilings, not requests. Each process asks for whatever its nf-core
  resource label specifies; only requests above these values are capped.

The retry strategy targets transient scheduler failures, not analysis errors.
Exit codes 104, 134, 137, 139, 143 and 247 cover crashes, out-of-memory kills
and SIGTERM from preemption. The presence of a `nonpre` partition on this
cluster implies that jobs in `main` may be preemptible, which makes the retry
worth having whether or not preemption is ever observed. Anything else uses
`finish` rather than `terminate`, so tasks already running are allowed to
complete and `-resume` has less to redo.

Container settings enable Singularity with `autoMounts`, which is required for
tasks to see `/scratch` from inside a container, and point the image cache at
scratch. Setting the cache directory in the config rather than relying on the
`NXF_SINGULARITY_CACHEDIR` environment variable keeps the behaviour with the
repository instead of in a shell profile.

---

## 2. Smoke test

Before committing the real samples, `workflow/run_taxprofiler_test.slurm` runs
the pipeline's own bundled test dataset with `-profile test,singularity`. It
takes about twenty minutes and separates two classes of failure that would
otherwise be tangled together: a broken cluster configuration, and a problem
with our data.

It also settles the one assumption in the design that could not be checked any
other way. The Nextflow head process runs as a SLURM job on a compute node, and
from there it must itself submit jobs. Most clusters permit this; not all do.
If Amarel did not, `executor = 'slurm'` would be unusable and the whole
configuration would need rethinking. The script tests for `sbatch` explicitly
before starting the pipeline, so that failure would surface in seconds rather
than twenty minutes in.

Result on 2026-08-31: 179 tasks succeeded, none failed, and the stderr log was
empty. The test profile exercises roughly fifteen different classifiers, so
every one of those was a separate container image pulled and executed
successfully. This is a stronger check of the container runtime than profiling
a single tool would have been.

The test script is kept in the repository rather than deleted. It is the
evidence that the configuration works, and it is the fastest way to re-verify
after any cluster change.

---

## 3. Production run

`workflow/run_taxprofiler.slurm` runs the 24 selected samples. The head process
requests one core, 8 GB and two days of wall clock: it submits work and waits,
so it needs time rather than resources. Two days is roughly ten times the
expected runtime, chosen so that a stalled run dies rather than holding a queue
slot for three days.

Analysis parameters and their basis:

| Parameter | Setting | Basis |
|---|---|---|
| run_metaphlan | on | The only classifier used; index pinned per docs/05 |
| perform_shortread_qc | on | fastp. Costs minutes and produces the MultiQC report |
| host removal | off | See below |
| perform_runmerging | off | Sample and run are one-to-one in this samplesheet |
| run_profile_standardisation | on | Produces the taxpasta merged table |
| metaphlan_save_samfiles | off (default) | SAM files would add tens of gigabytes with no downstream use |

Host removal is disabled on the evidence assembled in `docs/03`: across all
1,679 samples in the cohort, the ratio of trimmed to raw reads is unimodal at a
mean of 86.3% (min 67.1%, max 96.9%) and flat across sampling ages, with adult
maternal samples indistinguishable from newborns. Human DNA in infant stool
declines with age, so a substantial host fraction should have produced a
systematic gradient. It does not. The claim is not that the original authors
already removed host reads, which is not something this evidence can establish;
it is that the host fraction is too small to justify downloading a human
reference and building a bowtie2 index against it.

Result on 2026-08-31: 99 tasks succeeded, none failed, 38 minutes wall clock,
36 CPU hours. All four per-sample stages reported 24 of 24: FastQC before
trimming, fastp, FastQC after trimming, and MetaPhlAn.

---

## 4. Output structure

The results tree contains only the six directories corresponding to enabled
tools: `fastp`, `fastqc`, `metaphlan`, `multiqc`, `pipeline_info`, `taxpasta`.
Disabled classifiers produce no output at all, which is the visible confirmation
that the parameters took effect.

MetaPhlAn output is nested under the `db_name` from `data/databases.csv`, so the
path itself records which index version produced the results. Downstream code
must therefore read `db_name` from the database sheet or glob for it, rather
than hard-coding a directory name.

The merged table at `taxpasta/metaphlan_mpa_vJun23_202403.tsv` is the entry
point for downstream analysis: one file, all samples, 399 taxa. Two properties
of it matter:

- Abundances are stored as integers scaled by 1e8, not as percentages.
- Rows are cumulative across ranks. The row for taxonomy_id 2 (Bacteria) equals
  the sum of everything beneath it, so analysis must filter to a single rank
  before computing anything.

Column names take the form `{sample}_{run}_{db_name}.metaphlan_profile`. Since
sample and run hold the same accession in this samplesheet, the accession
appears twice and must be parsed out before use as a sample identifier.

Results are copied from scratch to `results/` in the repository directory after
each run (1.1 GB, ignored by git). Scratch was measured at 95% full on
2026-08-20 and is subject to purge policy; everything there is treated as
reconstructible, but re-running costs 38 minutes plus a 1.5 hour index rebuild
if the database is lost too.

---

## 5. Observation: fastp removes almost nothing

fastp retained 99.99997% of reads across all 24 samples, with adapter content
between 0.01% and 0.11%. Read counts before and after trimming are effectively
identical.

This is inconsistent with the 86.3% retention computed from the cohort metadata
in `docs/03`, which would predict roughly one read in seven being discarded.
The read length distributions explain the discrepancy:

| Measure | Value |
|---|---|
| Median read length, pre-trimming | 100 or 124, depending on sequencing batch |
| Mean read length, pre-trimming | 99.5, and 114.3 to 119.9 |

Mean below median means the length distribution already has a left tail.
Untrimmed Illumina reads are uniformly full length. The reads retrieved from ENA
therefore already carry the marks of quality trimming, and the published
raw-versus-trimmed counts describe processing that happened before deposition,
not processing that this pipeline performs.

The consequence is a scoping one rather than a correction: fastp here functions
as verification rather than as processing, and the 86.3% figure should not be
cited as describing this pipeline's own read handling. Keeping fastp enabled is
still worthwhile, since it is what established this fact and it supplies the
MultiQC report.

Two further observations from the same report, recorded because they are easy to
misread as faults:

- GC content ranges from 33% to 59% across samples. This tracks community
  composition rather than data quality: Bifidobacterium sits near 60% GC and
  Escherichia near 50%.
- FastQC and fastp report different duplication rates (median 21% versus 40%).
  They measure different things and are not expected to agree.

---

## Not yet done

- The taxonomic table has not been joined to birth mode and sampling age.
  Nothing in this document interprets the profiles biologically, and nothing
  should until that join is made and validated.
- Functional profiling with HUMAnN3 is not run. Its reference databases are not
  downloaded.
- Downstream analysis code does not exist yet.
