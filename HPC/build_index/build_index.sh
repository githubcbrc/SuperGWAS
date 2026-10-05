#!/usr/bin/env bash
# Build k-mer lookup index from merged kgwas results
# Direct binary invocation, no container.

set -euo pipefail

if [[ -z "${KGWAS_OUT:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

KGWAS_BATCH_DIR="${KGWAS_OUT}/${PHENO_NAME}"
OUTPUT="${INDEX_DIR}/${PHENO_NAME}"

[[ -d "$KGWAS_BATCH_DIR" ]] || { echo "ERROR: Directory not found: $KGWAS_BATCH_DIR" >&2; exit 1; }
TSV_COUNT=$(find "$KGWAS_BATCH_DIR" -maxdepth 1 -name '*.tsv' 2>/dev/null | wc -l)
[[ "$TSV_COUNT" -gt 0 ]] || { echo "ERROR: No .tsv files in $KGWAS_BATCH_DIR" >&2; exit 1; }
echo "Found $TSV_COUNT .tsv files in $KGWAS_BATCH_DIR"

mkdir -p "$OUTPUT"

echo "== build_index =="
echo "INPUT:  $KGWAS_BATCH_DIR"
echo "OUTPUT: $OUTPUT"

"${BIN_DIR}/build_index" "$KGWAS_BATCH_DIR" "$OUTPUT" --shards "${INDEX_SHARDS:-16}"

echo "Done. Index written to: $OUTPUT"
