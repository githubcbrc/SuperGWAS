#!/bin/bash
#SBATCH --job-name=kgwas_2d_merge
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

source "${REPO_ROOT}/config.sh"

OUT_PREFIX="${KGWAS_2D_OUT}/${PHENO_NAME}_2d"
NUM_2D_TASKS="${NUM_2D_TASKS:-1}"

echo "== kgwas_2d_merge =="
echo "Merging ${NUM_2D_TASKS} parts into ${OUT_PREFIX}.tsv"

# Concatenate all part files
cat "${OUT_PREFIX}".part_*.tsv > "${OUT_PREFIX}.tsv"

# Clean up parts
rm -f "${OUT_PREFIX}".part_*.tsv

echo "Done. $(wc -l < "${OUT_PREFIX}.tsv") total pairs"
echo "K-mers file: ${OUT_PREFIX}.kmers"
