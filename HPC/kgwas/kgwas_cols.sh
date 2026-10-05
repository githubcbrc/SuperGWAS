#!/usr/bin/env bash
# Worker script for kgwas --mode cols (direct binary, no container)
#
# Processes one batch identified by SLURM_ARRAY_TASK_ID.

set -euo pipefail

if [[ -z "${COL_MATRIX_DIR:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

for var in COL_MATRIX_DIR ACCESSIONS PHENOTYPE U_BIN KGWAS_OUT PHENO_NAME; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var not set" >&2
        exit 1
    fi
done

IDX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
BATCH="$(printf 'batch_%03d' "$IDX")"

BATCH_PATH="${COL_MATRIX_DIR}/${BATCH}"
OUT_DIR="${KGWAS_OUT}/${PHENO_NAME}"
OUT_TSV="${OUT_DIR}/${BATCH}.tsv"

LOG10_THRESH="${LOG10_THRESH:-6}"
CORR_THRESHOLD="${CORR_THRESHOLD:-0.2}"

for p in "${ACCESSIONS}" "${PHENOTYPE}" "${U_BIN}" "${BATCH_PATH}"; do
    if [[ ! -e "$p" ]]; then
        echo "[WARN] ${BATCH}: missing $p -- skipping."
        exit 0
    fi
done

mkdir -p "${OUT_DIR}"

NPROC="$(nproc)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-$NPROC}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export LOAD_THREADS="${LOAD_THREADS:-24}"
export BLAS_THREADS="${BLAS_THREADS:-$NPROC}"
export EMIT_THREADS="${EMIT_THREADS:-$NPROC}"
export CORR_THRESHOLD

echo "== kgwas --mode cols =="
echo "BATCH:          ${BATCH}"
echo "BATCH_PATH:     ${BATCH_PATH}"
echo "OUT_TSV:        ${OUT_TSV}"
echo "CORR_THRESHOLD: ${CORR_THRESHOLD}"
echo "LOG10_THRESH:   ${LOG10_THRESH}"

/usr/bin/time -v "${BIN_DIR}/kgwas" \
    "${ACCESSIONS}" \
    "${PHENOTYPE}" \
    "${BATCH_PATH}" \
    "${U_BIN}" \
    "${OUT_TSV}" \
    "${LOG10_THRESH}" \
    --mode cols
