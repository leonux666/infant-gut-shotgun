#!/bin/bash
#SBATCH --job-name=taxprofiler
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/xw347/infant-gut-shotgun/logs/taxprofiler-%j.out
#SBATCH --error=/home/xw347/infant-gut-shotgun/logs/taxprofiler-%j.err

# Taxonomic profiling of 24 infant gut shotgun metagenomes with MetaPhlAn 4.1.1.
#
# This job is the Nextflow head process only: it submits work to SLURM and waits.
# It needs wall time, not cores. Two days is roughly ten times the expected
# runtime; if it ever runs that long something is wrong and it should die rather
# than sit in the queue for three days.
#
# Analysis choices, each with its basis:
#   --run_metaphlan              The only classifier we profile with. Database is
#                                pinned to mpa_vJun23_CHOCOPhlAnSGB_202403.
#   --perform_shortread_qc       fastp. Costs minutes and produces the MultiQC
#                                report.
#   host removal                 NOT enabled. Read-retention across all 1,679
#                                samples is unimodal at 86.3% (min 67.1%, max
#                                96.9%) and flat across timepoints, with adult
#                                mothers indistinguishable from newborns. Human
#                                DNA is too small a fraction to justify building
#                                a bowtie2 index against the human genome.
#   --perform_runmerging         NOT enabled. Sample and run are one-to-one here.
#   --run_profile_standardisation
#                                Produces the taxpasta merged table, which is
#                                what downstream analysis reads.
#
# metaphlan_save_samfiles defaults to false, so no SAM files are written.

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
echo

echo "=== Inputs ==="
# Fail here rather than twenty minutes in: Nextflow validates db_path at startup
# but nothing checks that these two files exist before the run begins.
for f in data/samplesheet.csv data/databases.csv; do
    if [[ ! -s "$f" ]]; then
        echo "FATAL: $f is missing or empty"
        exit 1
    fi
    echo "$f: $(wc -l < "$f") lines"
done
echo

echo "=== Running taxprofiler ==="
nextflow run nf-core/taxprofiler \
    -r 2.0.1 \
    -profile singularity \
    -c "$PROJECT_DIR/config/amarel.config" \
    -w "$SCRATCH_DIR/work" \
    --input "$PROJECT_DIR/data/samplesheet.csv" \
    --databases "$PROJECT_DIR/data/databases.csv" \
    --outdir "$SCRATCH_DIR/results" \
    --perform_shortread_qc \
    --run_metaphlan \
    --run_profile_standardisation \
    -resume

echo
echo "=== Done ==="
echo "date: $(date)"
