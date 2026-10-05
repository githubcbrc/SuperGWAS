#!/usr/bin/env bash
# Worker script for kgwas --mode rows (direct binary, no container)
#
# Processes one chunk identified by SLURM_ARRAY_TASK_ID.

set -euo pipefail

if [[ -z "${ROW_MATRIX_DIR:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

for var in ROW_MATRIX_DIR ACCESSIONS PHENOTYPE U_BIN KGWAS_OUT PHENO_NAME; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var not set" >&2
        exit 1
    fi
done

IDX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"

CHUNK_FILE=$(ls -1 "${ROW_MATRIX_DIR}"/*.tsv* 2>/dev/null | sort -V | sed -n "$((IDX+1))p")

if [[ -z "$CHUNK_FILE" || ! -e "$CHUNK_FILE" ]]; then
    echo "[WARN] Chunk index ${IDX}: no file found -- skipping."
    exit 0
fi

CHUNK_NAME=$(basename "${CHUNK_FILE}" | sed 's/\.tsv\(\.gz\)\?$//')

OUT_DIR="${KGWAS_OUT}/${PHENO_NAME}"
OUT_TSV="${OUT_DIR}/${CHUNK_NAME}.tsv"
CORR_THRESHOLD="${CORR_THRESHOLD:-0.2}"

for p in "${ACCESSIONS}" "${PHENOTYPE}" "${U_BIN}" "${CHUNK_FILE}"; do
    if [[ ! -e "$p" ]]; then
        echo "[WARN] ${CHUNK_NAME}: missing $p -- skipping."
        exit 0
    fi
done

mkdir -p "${OUT_DIR}"

NPROC="$(nproc)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-$NPROC}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export CORR_THRESHOLD

echo "== kgwas --mode rows =="
echo "CHUNK_FILE:     ${CHUNK_FILE}"
echo "OUT_TSV:        ${OUT_TSV}"
echo "CORR_THRESHOLD: ${CORR_THRESHOLD}"

/usr/bin/time -v "${BIN_DIR}/kgwas" \
    "${ACCESSIONS}" \
    "${PHENOTYPE}" \
    "${CHUNK_FILE}" \
    "${U_BIN}" \
    "${OUT_TSV}" \
    --mode rows
