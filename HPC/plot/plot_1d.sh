#!/usr/bin/env bash
# Plot 1D Manhattan for a single reference (used as array job)
# Expects: SLURM_ARRAY_TASK_ID indexes into sorted ref directories
set -euo pipefail

if [[ -z "${KGWAS_OUT:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
module load python/3.11.0 2>/dev/null || true

KGWAS_BATCH_DIR="${KGWAS_OUT}/${PHENO_NAME}"
PLOT_1D_SCRIPT="${REPO_ROOT}/tools/plot_1d.py"

# Get the ref directory for this array task
mapfile -t REF_DIRS < <(find "${MAPPING_OUT}/${PHENO_NAME}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
IDX="${SLURM_ARRAY_TASK_ID:-0}"

if [[ "$IDX" -ge "${#REF_DIRS[@]}" ]]; then
    echo "Array index $IDX out of range (${#REF_DIRS[@]} refs)"
    exit 0
fi

REF_DIR_PATH="${REF_DIRS[$IDX]}"
REF_NAME=$(basename "$REF_DIR_PATH")
FILTERED="${REF_DIR_PATH}/${PHENO_NAME}.filtered.tsv"
PLOT_OUT_DIR="${REF_DIR_PATH}/plots"

if [[ ! -f "$FILTERED" ]]; then
    echo "SKIP: No filtered.tsv for ${REF_NAME}"
    exit 0
fi

mkdir -p "$PLOT_OUT_DIR"

REF_FASTA=$(find "${REF_DIR}" -maxdepth 1 \( -name "${REF_NAME}.fasta" -o -name "${REF_NAME}.fa" \) 2>/dev/null | head -1)

echo "Plotting 1D: ${REF_NAME}"
python3 "$PLOT_1D_SCRIPT" "$FILTERED" \
    --title "${PHENO_NAME} - ${REF_NAME}" \
    ${REF_FASTA:+--fasta "$REF_FASTA"} \
    --min-logp "${LOG10_MIN:-6}" \
    --out "$PLOT_OUT_DIR/${PHENO_NAME}.png" \
    --counts "$KGWAS_BATCH_DIR"/*.count 2>/dev/null

echo "Saved: $PLOT_OUT_DIR/${PHENO_NAME}.png"
