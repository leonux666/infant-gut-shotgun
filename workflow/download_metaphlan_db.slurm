#!/bin/bash
#SBATCH --job-name=mpa_db
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/xw347/infant-gut-shotgun/logs/mpa_db_%j.out
#SBATCH --error=/home/xw347/infant-gut-shotgun/logs/mpa_db_%j.err

source ~/.bashrc
conda activate nfcore

set -euo pipefail

DB_DIR=/scratch/xw347/nfcore/databases/metaphlan
INDEX=mpa_vJun23_CHOCOPhlAnSGB_202403
URL=https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/${INDEX}.tar
EXPECTED_BYTES=3316336640

cd "$DB_DIR"

echo "[1/4] downloading ${INDEX}.tar"
curl -sSL -o "${INDEX}.tar" "$URL"

ACTUAL_BYTES=$(stat -c %s "${INDEX}.tar")
echo "expected ${EXPECTED_BYTES} bytes, got ${ACTUAL_BYTES} bytes"
if [ "$ACTUAL_BYTES" -ne "$EXPECTED_BYTES" ]; then
    echo "ERROR: size mismatch, download incomplete"
    exit 1
fi

echo "[2/4] unpacking archive"
tar -xf "${INDEX}.tar"
ls -lh

echo "[3/4] decompressing SGB marker sequences"
bunzip2 -k "${INDEX}_SGB.fna.bz2"
echo "fasta headers: $(grep -c '^>' ${INDEX}_SGB.fna)"

echo "[4/4] building bowtie2 index"
bowtie2-build --threads 16 --large-index \
    "${INDEX}_SGB.fna" "${INDEX}"

echo "index files built:"
ls -lh ${INDEX}*.bt2l

echo "done"
