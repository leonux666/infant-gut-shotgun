# 01 — Environment Setup

Reproducible setup log for this project on the Rutgers Amarel HPC cluster
(Linux x86_64, Lmod module system, SLURM scheduler).

Every command below was run and verified. Commands are listed in the order
they should be executed on a fresh account.

---

## 1. Check for an existing conda installation

The cluster provides no conda module, so conda must be installed per user.
Check first — a previous installation may already exist.

```bash
# Look for a conda module (returns nothing on Amarel)
module avail 2>&1 | grep -i -E "conda|anaconda|miniforge|mamba"

# Check whether conda is already on PATH
conda --version
```

If `conda --version` prints a version, skip to step 3.

---

## 2. Install Miniforge (only if conda is absent)

Miniforge is preferred over Anaconda/Miniconda: it defaults to the
`conda-forge` channel, which — together with `bioconda` — hosts virtually all
bioinformatics packages, and it ships with `mamba` for faster dependency
resolution.

```bash
cd ~
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Verify the download is a complete shell script, not an error page
head -c 200 Miniforge3-Linux-x86_64.sh

# Verify payload integrity against the MD5 declared in the script header.
# The installer is a shell header concatenated with a compressed payload,
# so the declared checksum covers the payload only — not the whole file.
tail -n +$(( $(grep -anm 1 '^@@END_HEADER@@' Miniforge3-Linux-x86_64.sh | cut -d: -f1) + 1 )) Miniforge3-Linux-x86_64.sh | md5sum

# Install non-interactively into the home directory.
# Home is used rather than scratch: scratch is subject to purge policies,
# and the conda installation underpins every downstream analysis.
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3

# Write the shell initialisation block into ~/.bashrc
$HOME/miniforge3/bin/conda init bash

# Reload the shell configuration
source ~/.bashrc

# Clean up
rm ~/Miniforge3-Linux-x86_64.sh
```

---

## 3. Create the workflow environment

This environment holds only the workflow engine. Analysis software
(MetaPhlAn4, HUMAnN3, MetaBAT2, ...) is **not** installed here — nf-core
pipelines pull those as Singularity containers at runtime.
It also carries `pandas` and `openpyxl`, which the metadata scripts under `src/`
import directly. `pandas` arrives anyway as an indirect dependency of nf-core,
but an indirect dependency is not a declared one — it is pinned here so the
scripts do not silently break when nf-core changes its own requirements.

```bash
mamba create -n nfcore -c conda-forge -c bioconda nextflow nf-core pandas openpyxl -y
```

Channel order matters: `conda-forge` must precede `bioconda`, per bioconda's
own installation guidance.

Alternatively, recreate the environment from the file in this repository:

```bash
mamba env create -f environment.yml
```

---

## 4. Activate and verify

```bash
conda activate nfcore
nextflow -version
nf-core --version
```

Verified versions:

| Tool | Version |
| --- | --- |
| conda | 26.1.0 |
| mamba | 2.5.0 |
| Nextflow | 26.04.6 |
| nf-core/tools | 4.1.0 |
| pandas | 3.0.5 |
| openpyxl | 3.1.5 |

Nextflow requires a JVM. If it fails to start, load the cluster's Java module
first:

```bash
module load java
```

---

## 5. Export the environment file

```bash
conda env export --no-builds -n nfcore > environment.yml
```

`--no-builds` strips platform-specific build strings, which would otherwise
make the file unusable on any other machine.

The auto-generated export was replaced with a hand-written file listing only
direct dependencies with pinned versions. Two reasons: the full export runs
to ~200 lines of transitive dependencies that obscure intent, and it appends
a `prefix:` line containing the local installation path.

---

## Notes

- Analyses must **not** be run on Amarel login nodes. Package installation is
  I/O-bound and acceptable; anything compute-bound goes through SLURM.
- Containers are provided by Apptainer 1.5.3 at `/usr/bin/singularity` on the
  production login nodes. The `singularity/3.1.0` module must **not** be loaded:
  it is a non-functional leftover that shadows the working binary. See `docs/05`
  for the full diagnosis.
- Storage split: code and environments in `$HOME`; raw reads, intermediate
  files, and Nextflow work directories in `/scratch/$USER`.
