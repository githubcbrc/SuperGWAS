#!/usr/bin/env bash
# Submit convert_matrix step

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

export SCRATCH_ROOT="${SCRATCH_ROOT:-$RESULTS_ROOT}"
export SORTED_OUT="${SORTED_OUT:-${COL_MATRIX_DIR}}"

mkdir -p "${LOGS_DIR}/convert_matrix"

echo "Configuration:"
echo "  ROW_MATRIX_DIR: $ROW_MATRIX_DIR"
echo "  SORTED_OUT:     $SORTED_OUT"
echo "  SCRATCH_ROOT:   $SCRATCH_ROOT"
echo ""

"${SCRIPT_DIR}/convert_matrix.sbatch"
