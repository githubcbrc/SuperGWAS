#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------
# Submit Slurm array to process multiple batches in parallel.
# Direct binary invocation, no container.
# ----------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${REPO_ROOT}/config.sh"

# Settings
BATCH_SIZE="${BATCH_SIZE:-10}"
GLOB="${GLOB:-*.gz}"
SORTED_OUT="${SORTED_OUT:-${RESULTS_ROOT}/sorted}"
BATCH_TMP="${BATCH_TMP:-${RESULTS_ROOT}/batches}"

# 1) Prep batches
"${SCRIPT_DIR}/prep_batches.sh"

MANIFEST="${SORTED_OUT}/MANIFEST"
NUM_BATCHES="$(wc -l < "$MANIFEST")"
(( NUM_BATCHES > 0 )) || { echo "No batches listed in $MANIFEST" >&2; exit 1; }

ARRAY_RANGE="0-$((NUM_BATCHES-1))"
mkdir -p "${SORTED_OUT}/logs"

echo "[submit] ROW_MATRIX_DIR : $ROW_MATRIX_DIR"
echo "[submit] SORTED_OUT     : $SORTED_OUT"
echo "[submit] BATCH_TMP      : $BATCH_TMP"
echo "[submit] MANIFEST       : $MANIFEST ($NUM_BATCHES batches)"

TIME="${TIME:-04:00:00}"
CPUS="${CPUS:-48}"
MEM="${MEM:-256G}"

# 2) Submit array
sbatch \
  --job-name=convert_matrix \
  --array="$ARRAY_RANGE" \
  --time="$TIME" \
  --cpus-per-task="$CPUS" \
  --mem="$MEM" \
  --export=ALL \
  --output="${SORTED_OUT}/logs/%x_%A_%a.out" \
  --error="${SORTED_OUT}/logs/%x_%A_%a.err" \
  "${SCRIPT_DIR}/convert_matrix.sh"

echo "[submit] Submitted array: $ARRAY_RANGE"
echo "[submit] Logs: ${SORTED_OUT}/logs/"
