#!/bin/bash
#SBATCH --job-name=map_2d_pairs
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

set -euo pipefail

source "${REPO_ROOT}/config.sh"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

# Select reference (array job)
mapfile -t REF_FILES < <(find "${REF_DIR}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | sort)
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
REF_FASTA="${REF_FILES[$TASK_ID]}"
REF_NAME=$(basename "$REF_FASTA" | sed 's/\.\(fasta\|fa\)$//')

KMERS_FILE="${KGWAS_2D_OUT}/${PHENO_NAME}_2d.kmers"
PAIRS_FILE="${KGWAS_2D_OUT}/${PHENO_NAME}_2d.tsv"
OUT_DIR="${KGWAS_2D_OUT}/${REF_NAME}"
OUT_FILE="${OUT_DIR}/${PHENO_NAME}_2d.tsv"

mkdir -p "$OUT_DIR"

echo "== map_2d_pairs =="
echo "REF:     $REF_NAME"
echo "KMERS:   $KMERS_FILE"
echo "PAIRS:   $PAIRS_FILE"
echo "OUTPUT:  $OUT_FILE"

"${BIN_DIR}/map_2d_pairs" "$KMERS_FILE" "$PAIRS_FILE" "$REF_FASTA" "$K" "$OUT_FILE" \
    --threads "${SLURM_CPUS_PER_TASK:-16}"

echo "Done. Mapped pairs: $OUT_FILE"
