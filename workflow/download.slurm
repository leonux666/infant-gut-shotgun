#!/bin/bash
#SBATCH --job-name=fastq_dl
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err

# Fetch the 24 selected runs (48 files, ~29 GB) from the ENA mirror.
# The download is network-bound, so one core and little memory suffice.

# Load the environment before enabling strict mode: /etc/bashrc probes
# undefined variables, which set -u would treat as fatal.
source ~/.bashrc
conda activate nfcore

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

python src/download_fastq.py \
  --input data/selected_samples.tsv \
  --outdir data/fastq
