#!/usr/bin/env bash
# Plot 1D Manhattan + 2D 3D Manhattan for a single reference (array job)
# Expects: SLURM_ARRAY_TASK_ID indexes into sorted ref files
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

echo "== plot: ${REF_NAME} =="

# ===========================================================================
# 1D Plot
# ===========================================================================
KGWAS_BATCH_DIR="${KGWAS_OUT}/${PHENO_NAME}"
PLOT_1D_SCRIPT="${REPO_ROOT}/tools/plot_1d.py"
FILTERED="${MAPPING_OUT}/${PHENO_NAME}/${REF_NAME}/${PHENO_NAME}.filtered.tsv"
PLOT_OUT_DIR="${MAPPING_OUT}/${PHENO_NAME}/${REF_NAME}/plots"

if [[ -f "$FILTERED" ]]; then
    mkdir -p "$PLOT_OUT_DIR"
    echo "Plotting 1D: ${REF_NAME}"
    python3 "$PLOT_1D_SCRIPT" "$FILTERED" \
        --title "${PHENO_NAME} - ${REF_NAME}" \
        ${REF_FASTA:+--fasta "$REF_FASTA"} \
        --min-logp "${LOG10_MIN:-6}" \
        --out "$PLOT_OUT_DIR/${PHENO_NAME}.png" \
        --counts "$KGWAS_BATCH_DIR"/*.count 2>/dev/null || \
        echo "WARNING: 1D plot failed for ${REF_NAME}" >&2
    echo "Saved: $PLOT_OUT_DIR/${PHENO_NAME}.png"
else
    echo "SKIP 1D: No filtered.tsv for ${REF_NAME}"
fi

# ===========================================================================
# 2D Plot
# ===========================================================================
KGWAS_2D_TSV="${KGWAS_2D_OUT}/${REF_NAME}/${PHENO_NAME}_2d.tsv"
PLOT_2D_SCRIPT="${REPO_ROOT}/tools/plot_2d.py"
PLOT_2D_OUT="${KGWAS_2D_OUT}/${REF_NAME}/${PHENO_NAME}_2d.html"

if [[ -f "$KGWAS_2D_TSV" && -f "$PLOT_2D_SCRIPT" ]]; then
    echo "Plotting 2D: ${REF_NAME}"

    MARGINAL_FLAG=""
    [[ -f "$FILTERED" ]] && MARGINAL_FLAG="--marginal $FILTERED"

    EFFECT_FLAG=""
    [[ -n "${EFFECT_PERCENTILE:-}" ]] && EFFECT_FLAG="--effect-percentile $EFFECT_PERCENTILE"

    python3 "$PLOT_2D_SCRIPT" \
        --fasta "$REF_FASTA" \
        --tsv "$KGWAS_2D_TSV" \
        --title "${PHENO_NAME} - ${REF_NAME}" \
        $MARGINAL_FLAG \
        $EFFECT_FLAG \
        -o "$PLOT_2D_OUT" || \
        echo "WARNING: 2D plot failed for ${REF_NAME}" >&2
    echo "Saved: $PLOT_2D_OUT"
else
    echo "SKIP 2D: No 2D results for ${REF_NAME}"
fi

echo "== done: ${REF_NAME} =="
