#!/usr/bin/env bash
# Worker script for sample_kmers (direct binary, no container)
#
# Processes a batch of shards based on SLURM_ARRAY_TASK_ID.
# Uses central config file.

set -euo pipefail

# Load config if not already loaded
if [[ -z "${ROW_MATRIX_DIR:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

# Validate required variables
for var in ROW_MATRIX_DIR ACCESSIONS SAMPLED_OUT; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var not set" >&2
        exit 1
    fi
done

# Validate input
[[ -d "$ROW_MATRIX_DIR" ]] || { echo "ERROR: ROW_MATRIX_DIR not found: $ROW_MATRIX_DIR" >&2; exit 1; }
[[ -f "$ACCESSIONS" ]] || { echo "ERROR: ACCESSIONS not found: $ACCESSIONS" >&2; exit 1; }

# Create output directory
mkdir -p "$SAMPLED_OUT"

# Export environment for the binary
export RESERVOIR_SIZE="${RESERVOIR_SIZE:-10000000}"
export SHARDS_PER_JOB="${SHARDS_PER_JOB:-50}"
export SEED="${SEED:-42}"

echo "== sample_kmers =="
echo "ROW_MATRIX_DIR:  $ROW_MATRIX_DIR"
echo "ACCESSIONS:      $ACCESSIONS"
echo "SAMPLED_OUT:     $SAMPLED_OUT"
echo "RESERVOIR_SIZE:  $RESERVOIR_SIZE"
echo "SHARDS_PER_JOB:  $SHARDS_PER_JOB"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-0}"

"${BIN_DIR}/sample_kmers" "$ROW_MATRIX_DIR" "$ACCESSIONS" "$SAMPLED_OUT"

echo "Done."
