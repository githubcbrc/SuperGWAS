#!/bin/bash
#SBATCH --job-name=kgwas_2w
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G

set -euo pipefail

source "${REPO_ROOT}/config.sh"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

mkdir -p "$KGWAS_2W_OUT"

KGWAS_2W_TSV="${KGWAS_2W_OUT}/kgwas_2w.tsv"
EXTRA_FLAGS=""
[[ "${INTERACTION_ONLY:-false}" == "true" ]] && EXTRA_FLAGS="--interaction_only"

echo "== kgwas_2w =="
echo "U_2W:             $U_2W"
echo "TH_FILE:          $TH_FILE"
echo "TP_FILE:          $TP_FILE"
echo "HOST_MATRIX:      $HOST_MATRIX"
echo "PATHOGEN_MATRIX:  $PATHOGEN_MATRIX"
echo "PHENOTYPE_2W:     $PHENOTYPE_2W"
echo "OUT:              $KGWAS_2W_TSV"
[[ -n "$EXTRA_FLAGS" ]] && echo "FLAGS:            $EXTRA_FLAGS"

"${BIN_DIR}/kgwas_2w" $EXTRA_FLAGS \
    "$U_2W" "$TH_FILE" "$TP_FILE" \
    "$PHENOTYPE_2W" "$HOST_MATRIX" "$PATHOGEN_MATRIX" "$KGWAS_2W_TSV"

echo "Done. 2W results written to: $KGWAS_2W_TSV"
