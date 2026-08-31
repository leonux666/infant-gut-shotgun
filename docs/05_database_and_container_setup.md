# 05 — Reference database and container environment

Prepares the MetaPhlAn reference database and establishes a working container
runtime on the cluster. Both are prerequisites for running nf-core/taxprofiler;
neither is handled by the pipeline itself.

Two things are documented together here because they are causally linked: the
database had to be built by hand, with an explicit bowtie2 index step, partly
because no container runtime was available at the time to invoke MetaPhlAn's own
installer.

---

## 1. Container runtime

### The problem

`module load singularity` on Amarel provides Singularity 3.1.0, compiled
2019-02-25 against gcc-4.8. Running any image with it fails with
`FATAL: failed to resolve session directory`. The `/var/singularity` directory
that version requires does not exist on the cluster.

### Diagnosis

The decisive evidence was a prior success. The same QIIME2 image had run without
incident in March 2025 during 16S analysis; the exact command was still in shell
history. Re-running that command reproduced the failure.

This distinguishes between two hypotheses at almost no cost:

| Hypothesis | Prediction | Observed |
|---|---|---|
| Misconfiguration on our side | The old command would still work | It does not |
| The environment changed | The old command now fails | Confirmed |

Ruling out the first hypothesis redirected the search from "what did I set up
wrong" to "what changed on the cluster." The answer was in the SSH login banner,
not in any error message: the CentOS login nodes are being retired in favour of
Red Hat production nodes at `amarel.hpc.rutgers.edu`.

`amarel4` is one of the decommissioned CentOS nodes. `/var/singularity`
disappeared when it stopped being maintained. The `singularity/3.1.0` module is
a leftover pointing at a runtime that no longer has its supporting directories.

### Resolution

| | Old node (amarel4, CentOS) | Production nodes (amarel.hpc.rutgers.edu) |
|---|---|---|
| module load singularity | 3.1.0, non-functional | Same stale module, still listed |
| /var/singularity | absent | absent |
| /usr/bin/singularity | — | Apptainer 1.5.3, system-provided |

Two rules follow:

1. Connect to `amarel.hpc.rutgers.edu`, never `amarel4`.
2. Do not run `module load singularity`. The system already provides Apptainer
   1.5.3 under the compatibility name `singularity`; loading the module shadows
   it with the broken 3.1.0 binary.

nf-core detects Apptainer natively, so `-profile singularity` works unchanged.

Verified end to end on 2026-08-31: a full `-profile test,singularity` run of
nf-core/taxprofiler 2.0.1 completed 179 tasks across roughly fifteen distinct
tool containers, with zero failures and an empty stderr log.

---

## 2. MetaPhlAn database

### Version selection

MetaPhlAn database names carry two independent version components: a marker set
version and a build date, as in `mpa_v{markers}_CHOCOPhlAnSGB_{build}`.

The pipeline's MetaPhlAn version is 4.1.1, read from the container declaration
in `modules/nf-core/metaphlan/metaphlan/main.nf`. The pipeline's environment
file pins the tool version but says nothing about which index to pair with it.

`mpa_vJun23_CHOCOPhlAnSGB_202403` was selected on three independent grounds:

1. The MetaPhlAn 4.1.1 release notes name it as the paired index.
2. Independent users report running 4.1.1 against it successfully.
3. The changelog records that vJun23_202403 contains the same SGBs as
   vJun23_202307, differing only in corrected NCBI taxonomic assignments, which
   improves relative abundance estimates at higher ranks. Same content, with a
   bug fixed.

The vJan25 and vJan26 indexes target MetaPhlAn 4.2 and later; MetaPhlAn refuses
to run against a mismatched index.

The index must be pinned explicitly. Left to itself, MetaPhlAn queries
`mpa_latest` at runtime and will interrupt a run to ask whether to download a
newer index. An index that changes underneath a pipeline is not reproducible.

### Why not `metaphlan --install`

No callable MetaPhlAn binary exists outside the pipeline's containers: the
`nfcore` conda environment deliberately omits analysis tools, since nf-core
manages their versions itself. The installer performs download, unpack,
decompress, and index build. Writing those four stages as an explicit SLURM
script makes each one observable, verifiable, and documentable, which is worth
more here than the convenience of a single command.

### Build

`workflow/download_metaphlan_db.slurm` (16 cores, 64 GB, 24 h wall clock) runs
four stages, each printing a checkable number. Submit it with `sbatch`.

Measured facts:

| Item | Value |
|---|---|
| Source host | cmprod1.cibio.unitn.it, biobakery4 database directory |
| Protocol | HTTPS required; HTTP returns a 301 that curl ignores without -L |
| Tar size | 3,316,336,640 bytes, verified as a constant in the script |
| Tar contents | .pkl (94 M), _SGB.fna.bz2 (2.9 G), _VSG.fna.bz2 (174 M), _VINFO.csv (44 K) |
| Prebuilt index | Not included — inferable from tar size alone, since an index would push it past 20 GB |
| SGB sequences | 7,339,971 |
| Index output | Six .bt2l files, about 20 GB |
| bowtie2-build runtime | About 1.5 h on 16 cores |

Two details that are easy to get wrong:

- `bowtie2-build` requires `--large-index` to emit `.bt2l`. The combined
  sequence length is far beyond the 4 Gbp threshold for the 32-bit format.
- The index basename must match the `.pkl` basename, not the input FASTA name,
  which carries an extra `_SGB` suffix. MetaPhlAn constructs both filenames from
  the `--index` value, and a mismatch surfaces as a missing-database error that
  gives no hint that the cause is a filename.

Only the SGB index is built. The VSG file holds viral sequences used only under
`--profile_vsc`, which this project does not enable.

---

## 3. Database sheet

taxprofiler takes a second CSV alongside the samplesheet. The authoritative
column definition is `assets/schema_database.json` in the pipeline repository,
not the web documentation, which is inconsistent on this point. Read it from a
shallow clone of the tagged release rather than from the docs site.

Five columns, three required:

| Column | Required | Constraint |
|---|---|---|
| tool | yes | Enum of 14 values. MetaPhlAn is `metaphlan`, not `metaphlan4` |
| db_name | yes | No whitespace permitted |
| db_params | no | Extra arguments, unquoted. Empty is valid |
| db_type | no | short / long / short;long |
| db_path | yes | Declared with "exists": true, so it is validated at startup |

A uniqueEntries constraint applies to the tool and db_name pair, so a second
database can later be added under a different tool without renaming.

`data/databases.csv` is not tracked in git: it is short and hand-written, but
db_path is an absolute path valid only on this cluster. It holds a single row
naming the metaphlan tool, the db_name `mpa_vJun23_202403`, an empty db_params,
db_type short, and the absolute path to the database directory on scratch.

- db_params is empty, but both delimiting commas must still be present.
- db_type is set to short explicitly rather than accepting the short;long
  default, since all reads here are Illumina short reads.
- db_name becomes a directory name in the results tree, so encoding the index
  version in it makes the output self-describing.
- Check for stray carriage returns with `cat -A` before use: a single ^M from a
  copy-paste is enough to break parsing.

---

## Storage

| Path | Contents | Size |
|---|---|---|
| /scratch/xw347/nfcore/databases/metaphlan/ | Database and bowtie2 index | about 37 GB |
| /scratch/xw347/nfcore/singularity/ | Container image cache | grows per run |

Neither is tracked in git. Both are reconstructible: the database from the SLURM
script above, the images from the pipeline itself.

About 14 GB inside the database directory is recoverable once the index is
built, namely the source tar and the decompressed SGB FASTA. The compressed
FASTA is worth keeping, since it allows an index rebuild without re-downloading.

---

## Not yet done

- HUMAnN3 requires its own reference databases (ChocoPhlAn and UniRef90, tens of
  gigabytes). Not downloaded. Scratch was measured at 95% full on 2026-08-20,
  so this needs a capacity check before it is attempted.
