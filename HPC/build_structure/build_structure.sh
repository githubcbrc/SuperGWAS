#!/usr/bin/env bash
# Build structure input and compute population structure (U.bin)
# Direct binary invocation, no container.
#
# Steps:
#   1. Concat sample partitions (if distributed sampling was used)
#   2. Build structure input matrix (X_sample.tsv)
#   3. Compute population structure (U.bin)

set -euo pipefail

# Load config if not already loaded
if [[ -z "${ACCESSIONS:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

# Validate required variables
for var in SAMPLED_OUT SAMPLED_KMERS ACCESSIONS PHENOTYPE STRUCTURE_OUT U_BIN RANK; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var not set" >&2
        exit 1
    fi
done

for p in "$ACCESSIONS" "$PHENOTYPE"; do
    [[ -e "$p" ]] || { echo "ERROR: File not found: $p" >&2; exit 1; }
done

mkdir -p "$STRUCTURE_OUT"

# ===========================================================================
# Concat partitions if SAMPLED_KMERS doesn't exist
# ===========================================================================
if [[ ! -f "$SAMPLED_KMERS" ]]; then
    echo "== Concatenating sample partitions =="

    PARTS=($(ls -1 "${SAMPLED_OUT}"/part_*.tsv.gz 2>/dev/null | sort -V))

    if [[ ${#PARTS[@]} -eq 0 ]]; then
        echo "ERROR: SAMPLED_KMERS not found and no part_*.tsv.gz in $SAMPLED_OUT" >&2
        exit 1
    fi

    echo "Found ${#PARTS[@]} partition files"
    echo "Concatenating to ${SAMPLED_KMERS}..."

    TEMP_OUT="${SAMPLED_KMERS}.tmp.$$"
    cat "${PARTS[@]}" > "${TEMP_OUT}"
    mv "${TEMP_OUT}" "${SAMPLED_KMERS}"

    TOTAL_SEEN=0
    TOTAL_KEPT=0
    for stats in "${SAMPLED_OUT}"/part_*.stats; do
        if [[ -f "$stats" ]]; then
            seen=$(grep "^total_seen" "$stats" 2>/dev/null | cut -f2 || echo 0)
            kept=$(grep "^kept" "$stats" 2>/dev/null | cut -f2 || echo 0)
            TOTAL_SEEN=$((TOTAL_SEEN + seen))
            TOTAL_KEPT=$((TOTAL_KEPT + kept))
        fi
    done

    if [[ $TOTAL_SEEN -gt 0 ]]; then
        echo "Total k-mers seen: ${TOTAL_SEEN}"
        echo "Total k-mers kept: ${TOTAL_KEPT}"
    fi

    echo "Done concatenating: ${SAMPLED_KMERS}"
    echo ""
fi

[[ -f "$SAMPLED_KMERS" ]] || { echo "ERROR: SAMPLED_KMERS not found: $SAMPLED_KMERS" >&2; exit 1; }

X_SAMPLE="${STRUCTURE_OUT}/X_sample.tsv"

echo "== filter_kmers =="
echo "SAMPLED_KMERS:  $SAMPLED_KMERS"
echo "ACCESSIONS:     $ACCESSIONS"
echo "PHENOTYPE:      $PHENOTYPE"
echo "X_SAMPLE:       $X_SAMPLE"
echo "SUBSAMPLE_SIZE: $SUBSAMPLE_SIZE"
echo "MAF_THRESHOLD:  $MAF_THRESHOLD"
echo "SEED:           $SEED"

export SUBSAMPLE_SIZE MAF_THRESHOLD SEED

"${BIN_DIR}/filter_kmers" "$SAMPLED_KMERS" "$ACCESSIONS" "$PHENOTYPE" "$X_SAMPLE"

echo "Done. X_sample.tsv written to: $X_SAMPLE"

# ===========================================================================
# Step 3: Compute structure (U.bin + U2D.bin)
# ===========================================================================
echo ""
echo "== compute_structure_all =="
echo "X_SAMPLE:       $X_SAMPLE"
echo "PHENOTYPE:      $PHENOTYPE"
echo "RANK:           $RANK"
echo "U_BIN:          $U_BIN"
echo "U2D_BIN:        ${U2D_BIN:-${STRUCTURE_2D_OUT}/U2D.bin}"

mkdir -p "${STRUCTURE_2D_OUT:-${STRUCTURE_OUT}}"

"${BIN_DIR}/compute_structure_all" "$X_SAMPLE" "$PHENOTYPE" "$RANK" "$STRUCTURE_OUT" 1

# compute_structure_all writes U.bin and U2D.bin to STRUCTURE_OUT
# Copy U2D.bin to STRUCTURE_2D_OUT if it's a different directory
if [[ "${STRUCTURE_2D_OUT}" != "${STRUCTURE_OUT}" ]]; then
    cp "${STRUCTURE_OUT}/U2D.bin" "${STRUCTURE_2D_OUT}/U2D.bin"
fi

echo "Done. U.bin + U2D.bin written to: $STRUCTURE_OUT"
