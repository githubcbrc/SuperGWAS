#!/usr/bin/env bash
# =============================================================================
# Filter 1D kgwas batch results by p-value and merge into a single file.
#
# Usage:
#   ./filter_merge.sh [MAX_PVAL]
#   ./filter_merge.sh 1e-15      # keep only p < 1e-15 (default)
#   ./filter_merge.sh 1e-6       # keep only p < 1e-6
#   ./filter_merge.sh --dry-run  # show what would be submitted
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${REPO_ROOT}/config.sh"

MAX_PVAL="1e-15"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)   DRY_RUN=true; shift ;;
        -h|--help)   head -12 "$0" | tail -9; exit 0 ;;
        -*)          echo "Unknown option: $1" >&2; exit 1 ;;
        *)           MAX_PVAL="$1"; shift ;;
    esac
done

KGWAS_DIR="${KGWAS_OUT}/${PHENO_NAME}"
MERGED="${KGWAS_DIR}/${PHENO_NAME}.all.tsv"
FILTERED_DIR="${KGWAS_DIR}/filtered_tmp"

NUM_BATCHES=$(ls -1 "${KGWAS_DIR}"/batch_*.tsv 2>/dev/null | wc -l)
if [[ "$NUM_BATCHES" -eq 0 ]]; then
    echo "ERROR: No batch_*.tsv files in ${KGWAS_DIR}" >&2
    exit 1
fi

mkdir -p "${LOGS_DIR}/filter_merge"

echo "=============================================="
echo "Filter & Merge 1D Results"
echo "=============================================="
echo "Input dir:   ${KGWAS_DIR}"
echo "Batches:     ${NUM_BATCHES}"
echo "Max p-value: ${MAX_PVAL}"
echo "Output:      ${MERGED}"
echo "Dry run:     ${DRY_RUN}"
echo "=============================================="
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY-RUN] Would submit single job with ${NUM_BATCHES}-way parallel filter + merge"
    exit 0
fi

CPUS=48

JOB_ID=$(sbatch --parsable \
    --job-name=filter_merge \
    --time=00:30:00 \
    --cpus-per-task=${CPUS} \
    --mem=16G \
    -o "${LOGS_DIR}/filter_merge/%x.%j.out" \
    -e "${LOGS_DIR}/filter_merge/%x.%j.err" \
    --wrap "
set -euo pipefail
FILTERED_DIR='${FILTERED_DIR}'
mkdir -p \"\${FILTERED_DIR}\"

ls -1 ${KGWAS_DIR}/batch_*.tsv \
    | xargs -P ${CPUS} -I{} bash -c '
        base=\$(basename \"{}\")
        awk -F\"\t\" \"\\\$4 < ${MAX_PVAL}\" \"{}\" > \"'\${FILTERED_DIR}'/\${base%.tsv}.filtered.tsv\"
    '

echo \"Filtered ${NUM_BATCHES} batches (p < ${MAX_PVAL})\"

cat \${FILTERED_DIR}/batch_*.filtered.tsv > ${MERGED}
echo \"Merged to ${MERGED}: \$(wc -l < ${MERGED}) total rows\"
rm -rf \"\${FILTERED_DIR}\"
")

echo "[SUBMIT] filter_merge: ${JOB_ID}"
echo ""
echo "After completion, run the 2D pipeline:"
echo "  ./scripts/run_pipeline_2d.sh"
