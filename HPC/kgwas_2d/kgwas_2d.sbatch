#!/bin/bash
#SBATCH --job-name=kgwas_2d
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G

set -euo pipefail

source "${REPO_ROOT}/config.sh"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

INDEX_2D_DIR="${INDEX_2D_DIR:-${INDEX_DIR}/${PHENO_NAME}_2d}"
OUT_PREFIX="${KGWAS_2D_OUT}/${PHENO_NAME}_2d"
NUM_2D_TASKS="${NUM_2D_TASKS:-1}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

mkdir -p "$KGWAS_2D_OUT"

echo "== kgwas_2d (task ${TASK_ID}/${NUM_2D_TASKS}) =="
echo "INDEX_2D:      $INDEX_2D_DIR"
echo "U2D_BIN:       $U2D_BIN"
echo "LOG10_THRESH:  $LOG10_THRESH_2D"
echo "OUT_PREFIX:    $OUT_PREFIX"

"${BIN_DIR}/kgwas_2d" "$ACCESSIONS" "$PHENOTYPE" "$U2D_BIN" \
    "$INDEX_2D_DIR" "$OUT_PREFIX" "$LOG10_THRESH_2D" \
    --task-id "$TASK_ID" --num-tasks "$NUM_2D_TASKS"

echo "Done. Task ${TASK_ID} complete."
