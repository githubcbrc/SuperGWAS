#!/usr/bin/env bash
# Plot 2D 3D Manhattan for a single reference (array job support)
set -euo pipefail

if [[ -z "${KGWAS_OUT:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
module load python/3.11.0 2>/dev/null || true

# Select reference (array job support)
mapfile -t REF_FILES < <(find "${REF_DIR}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | sort)
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
REF_FASTA="${REF_FILES[$TASK_ID]}"
REF_NAME=$(basename "$REF_FASTA" | sed 's/\.\(fasta\|fa\)$//')

KGWAS_2D_TSV="${KGWAS_2D_OUT}/${REF_NAME}/${PHENO_NAME}_2d.tsv"
PLOT_2D_SCRIPT="${REPO_ROOT}/tools/plot_2d.py"
PLOT_2D_OUT="${KGWAS_2D_OUT}/${REF_NAME}/${PHENO_NAME}_2d.html"

if [[ ! -f "$KGWAS_2D_TSV" || ! -f "$PLOT_2D_SCRIPT" ]]; then
    echo "SKIP: no 2D results for ${REF_NAME}"
    exit 0
fi

echo "== plot_2d: ${REF_NAME} =="
echo "INPUT:  $KGWAS_2D_TSV"
echo "FASTA:  $REF_FASTA"
echo "OUTPUT: $PLOT_2D_OUT"

# Find 1D marginal for floor projection (same ref)
MARGINAL_1D="${MAPPING_OUT}/${PHENO_NAME}/${REF_NAME}/${PHENO_NAME}.filtered.tsv"
MARGINAL_FLAG=""
[[ -f "$MARGINAL_1D" ]] && MARGINAL_FLAG="--marginal $MARGINAL_1D"

EFFECT_FLAG=""
[[ -n "${EFFECT_PERCENTILE:-}" ]] && EFFECT_FLAG="--effect-percentile $EFFECT_PERCENTILE"

python3 "$PLOT_2D_SCRIPT" \
    --fasta "$REF_FASTA" \
    --tsv "$KGWAS_2D_TSV" \
    --title "${PHENO_NAME} - ${REF_NAME}" \
    $MARGINAL_FLAG \
    $EFFECT_FLAG \
    -o "$PLOT_2D_OUT"

echo "Saved: $PLOT_2D_OUT"
