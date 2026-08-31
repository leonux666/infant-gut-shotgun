#!/bin/bash
#SBATCH --job-name=taxprof-test
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/home/xw347/infant-gut-shotgun/logs/taxprofiler-test-%j.out
#SBATCH --error=/home/xw347/infant-gut-shotgun/logs/taxprofiler-test-%j.err

# Smoke test for the Amarel cluster configuration.
#
# This runs nf-core/taxprofiler against its own bundled test dataset (a few MB),
# not against our samples. Its only purpose is to prove four things at once:
#   1. The Nextflow head process can submit sbatch jobs from a compute node.
#   2. -profile singularity resolves to Apptainer 1.5.3 and pulls images.
#   3. config/amarel.config parses and its executor/queue settings take effect.
#   4. Images land in NXF_SINGULARITY_CACHEDIR on scratch, not in $HOME.
#
# Wall clock: expect 15-40 minutes, most of it pulling container images.
#
# This job is the Nextflow head process only. It uses almost no CPU; it submits
# work to the scheduler and waits. What it needs is wall time, not cores.

source ~/.bashrc

# Strict mode goes AFTER sourcing .bashrc: the system /etc/bashrc probes
# undefined variables, and set -u would kill the job before it starts.
set -euo pipefail

conda activate nfcore

PROJECT_DIR="/home/xw347/infant-gut-shotgun"
SCRATCH_DIR="/scratch/xw347/nfcore"

cd "$PROJECT_DIR"

echo "=== Environment ==="
echo "host:        $(hostname)"
echo "date:        $(date)"
echo "nextflow:    $(nextflow -v)"
echo "singularity: $(singularity --version)"
echo "NXF_SINGULARITY_CACHEDIR: ${NXF_SINGULARITY_CACHEDIR:-<unset>}"
echo

echo "=== Can this compute node submit jobs? ==="
# If sbatch is unavailable here, the whole slurm-executor design fails, and it
# is far cheaper to find out on this line than 20 minutes into a pipeline run.
command -v sbatch || { echo "FATAL: sbatch not available on compute nodes"; exit 1; }
echo

echo "=== Running taxprofiler test profile ==="
nextflow run nf-core/taxprofiler \
    -r 2.0.1 \
    -profile test,singularity \
    -c "$PROJECT_DIR/config/amarel.config" \
    -w "$SCRATCH_DIR/work-test" \
    --outdir "$SCRATCH_DIR/results-test" \
    -resume

echo
echo "=== Done ==="
echo "date: $(date)"
