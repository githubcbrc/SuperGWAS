#!/bin/bash
#SBATCH --job-name=compute_structure_2w
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

source "${REPO_ROOT}/config.sh"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

mkdir -p "$STRUCTURE_2W_OUT"

echo "== compute_structure_2w =="
echo "HOST_MATRIX:      $HOST_MATRIX"
echo "PATHOGEN_MATRIX:  $PATHOGEN_MATRIX"
echo "PHENOTYPE_2W:     $PHENOTYPE_2W"
echo "RANK_HOST:        $RANK_HOST"
echo "RANK_PATHOGEN:    $RANK_PATHOGEN"
echo "OUT:              $STRUCTURE_2W_OUT"

"${BIN_DIR}/compute_structure_2w" \
    "$PHENOTYPE_2W" "$HOST_MATRIX" "$PATHOGEN_MATRIX" \
    "$RANK_HOST" "$RANK_PATHOGEN" "${STRUCTURE_2W_OUT}/"

echo "Done. U0.txt written to: $U_2W"
