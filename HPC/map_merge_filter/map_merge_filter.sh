#!/usr/bin/env bash
# Worker script for map_merge_filter (direct binary, no container)
#
# Runs mapping, merging, and filtering for one reference genome.
# Uses SLURM_ARRAY_TASK_ID to select reference from REF_DIR.

set -euo pipefail

if [[ -z "${KGWAS_OUT:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

for var in KGWAS_OUT PHENO_NAME INDEX_DIR REF_DIR MAPPING_OUT K; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var not set" >&2
        exit 1
    fi
done

[[ -d "$REF_DIR" ]] || { echo "ERROR: REF_DIR not found: $REF_DIR" >&2; exit 1; }

mapfile -t REF_FILES < <(find "${REF_DIR}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | sort)
NUM_REFS=${#REF_FILES[@]}

if [[ "$NUM_REFS" -eq 0 ]]; then
    echo "ERROR: No .fasta or .fa files found in $REF_DIR" >&2
    exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

if [[ "$TASK_ID" -ge "$NUM_REFS" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID ($TASK_ID) >= number of refs ($NUM_REFS)" >&2
    exit 1
fi

REF_FASTA="${REF_FILES[$TASK_ID]}"
REF_NAME=$(basename "$REF_FASTA" | sed 's/\.\(fasta\|fa\)$//')

INDEX_PATH="${INDEX_DIR}/${PHENO_NAME}"
OUTPUT_DIR="${MAPPING_OUT}/${PHENO_NAME}/${REF_NAME}"
OUTPUT_FILE="${OUTPUT_DIR}/${PHENO_NAME}.filtered.tsv"

[[ -d "$INDEX_PATH" ]] || { echo "ERROR: Index not found: $INDEX_PATH" >&2; exit 1; }
[[ -f "$REF_FASTA" ]] || { echo "ERROR: Reference not found: $REF_FASTA" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"

NPROC="$(nproc)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-$NPROC}"

# Optional 2D output
EMIT_2D_ARGS=""
if [[ -n "${KGWAS_2D_OUT:-}" ]]; then
    EMIT_2D_FILE="${KGWAS_2D_OUT}/${PHENO_NAME}_2d_input.tsv"
    EMIT_2D_ARGS="--emit-2d ${EMIT_2D_FILE} --min-logp-2d ${LOG10_THRESH_2D_INPUT:-15}"
    mkdir -p "${KGWAS_2D_OUT}"
fi

echo "== map_merge_filter =="
echo "REF_NAME:     $REF_NAME"
echo "REF_FASTA:    $REF_FASTA"
echo "INDEX_PATH:   $INDEX_PATH"
echo "K:            $K"
echo "OUTPUT:       $OUTPUT_FILE"
echo "THREADS:      $OMP_NUM_THREADS"
echo "LOG10_MIN:    ${LOG10_MIN:-6}"
echo "WINDOW_SIZE:  ${WINDOW_SIZE:-100}"
echo "MIN_DENSITY:  ${MIN_DENSITY:-0.95}"
echo "CHUNK_SIZE:   ${CHUNK_SIZE:-10000000}"
if [[ -n "$EMIT_2D_ARGS" ]]; then
    echo "EMIT_2D:      $EMIT_2D_FILE"
    echo "MIN_LOGP_2D:  ${LOG10_THRESH_2D_INPUT:-15}"
fi

"${BIN_DIR}/map_merge_filter" \
    "$INDEX_PATH" \
    "$REF_FASTA" \
    "$K" \
    "$OUTPUT_FILE" \
    --threads "$OMP_NUM_THREADS" \
    --min-logp "${LOG10_MIN:-6}" \
    --window-size "${WINDOW_SIZE:-100}" \
    --min-density "${MIN_DENSITY:-0.95}" \
    --chunk-size "${CHUNK_SIZE:-10000000}" \
    $EMIT_2D_ARGS

echo "Done mapping. Output: $OUTPUT_FILE"
