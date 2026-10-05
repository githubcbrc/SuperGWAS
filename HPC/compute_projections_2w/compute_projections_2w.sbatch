#!/bin/bash
#SBATCH --job-name=compute_projections_2w
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

source "${REPO_ROOT}/config.sh"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

echo "== compute_projections_2w =="
echo "U_2W:             $U_2W"
echo "HOST_MATRIX:      $HOST_MATRIX"
echo "PATHOGEN_MATRIX:  $PATHOGEN_MATRIX"
echo "PHENOTYPE_2W:     $PHENOTYPE_2W"

"${BIN_DIR}/compute_projections_2w" \
    "$U_2W" "$PHENOTYPE_2W" "$HOST_MATRIX" "$PATHOGEN_MATRIX" "${STRUCTURE_2W_OUT}/"

echo "Done. Th.txt and Tp.txt written to: $STRUCTURE_2W_OUT"
