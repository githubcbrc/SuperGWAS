#!/usr/bin/env bash
# Build small 2D index from .2d.tsv files (k-mers meeting 2D p-value threshold)
set -euo pipefail

if [[ -z "${KGWAS_OUT:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/../../config.sh"
fi

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

KGWAS_2D_TSV_DIR="${KGWAS_OUT}/${PHENO_NAME}"
OUTPUT="${INDEX_2D_DIR:-${INDEX_DIR}/${PHENO_NAME}_2d}"

# Validate .2d.tsv files exist
TSV_2D_COUNT=$(find "$KGWAS_2D_TSV_DIR" -maxdepth 1 -name '*.2d.tsv' 2>/dev/null | wc -l)
[[ "$TSV_2D_COUNT" -gt 0 ]] || { echo "ERROR: No .2d.tsv files in $KGWAS_2D_TSV_DIR" >&2; exit 1; }
echo "Found $TSV_2D_COUNT .2d.tsv files in $KGWAS_2D_TSV_DIR"

mkdir -p "$OUTPUT"

# build_index expects .tsv extension — symlink .2d.tsv files
TMPDIR_2D=$(mktemp -d)
trap "rm -rf $TMPDIR_2D" EXIT
for f in "${KGWAS_2D_TSV_DIR}"/*.2d.tsv; do
    base=$(basename "$f" .2d.tsv)
    ln -s "$f" "${TMPDIR_2D}/${base}.tsv"
done

echo "== build_index_2d =="
echo "INPUT:  $KGWAS_2D_TSV_DIR (*.2d.tsv)"
echo "OUTPUT: $OUTPUT"

"${BIN_DIR}/build_index" "$TMPDIR_2D" "$OUTPUT" --shards 4

echo "Done. 2D index written to: $OUTPUT"
