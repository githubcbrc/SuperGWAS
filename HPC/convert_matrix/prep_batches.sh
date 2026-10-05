#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------
# Prepare batch staging dirs for convert_matrix
# -----------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${ROW_MATRIX_DIR:-}" ]]; then
    source "${SCRIPT_DIR}/../../config.sh"
fi

BATCH_SIZE="${BATCH_SIZE:-10}"
GLOB="${GLOB:-*.gz}"
SORTED_OUT="${SORTED_OUT:-${RESULTS_ROOT}/sorted}"
BATCH_TMP="${BATCH_TMP:-${RESULTS_ROOT}/batches}"

[[ -d "$ROW_MATRIX_DIR" ]] || { echo "ROW_MATRIX_DIR not found: $ROW_MATRIX_DIR" >&2; exit 1; }
mkdir -p "$BATCH_TMP" "$SORTED_OUT"

mapfile -d '' FILES < <(find "$ROW_MATRIX_DIR" -maxdepth 1 -type f -name "$GLOB" -print0 | sort -zV)
COUNT=${#FILES[@]}
(( COUNT > 0 )) || { echo "No files matched '$GLOB' in $ROW_MATRIX_DIR" >&2; exit 1; }

MANIFEST="${SORTED_OUT}/MANIFEST"
: > "$MANIFEST"

idx=0
batch=0
while (( idx < COUNT )); do
  bname=$(printf "batch_%03d" "$batch")
  bdir="$BATCH_TMP/$bname"
  outdir="$SORTED_OUT/$bname"

  rm -rf "$bdir"
  mkdir -p "$bdir" "$outdir"

  for f in "${FILES[@]:idx:BATCH_SIZE}"; do
    base="$(basename "$f")"
    ln -f "$f" "$bdir/$base" 2>/dev/null || ln -sfn "$f" "$bdir/$base"
  done

  echo "$bname" >> "$MANIFEST"
  (( batch++, idx+=BATCH_SIZE ))
done

printf '[prep] ROW_MATRIX_DIR  : %s\n' "$ROW_MATRIX_DIR"
printf '[prep] BATCH_TMP   : %s\n' "$BATCH_TMP"
printf '[prep] SORTED_OUT  : %s\n' "$SORTED_OUT"
printf '[prep] MANIFEST    : %s (%d batches)\n' "$MANIFEST" "$batch"
