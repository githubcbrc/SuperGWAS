#!/usr/bin/env bash
# Usage: merge.sh PHENO IPATH OPATH
# Concatenate "$IPATH"/*.tsv into "$OPATH/$PHENO.all.tsv"
# Also merges .count files for Bonferroni correction

set -euo pipefail

if (( $# != 3 )); then
  echo "Usage: $0 PHENO IPATH OPATH" >&2
  exit 2
fi

PHENO=$1
IPATH=$2
OPATH=$3

mkdir -p -- "$OPATH"
out="$OPATH/$PHENO.all.tsv"

cat "$IPATH"/*.tsv > "$out"

echo "Wrote: $out"

# Merge .count files (for Bonferroni correction)
countout="$OPATH/$PHENO.all.counts"
if ls "$IPATH"/*.count >/dev/null 2>&1; then
  cat "$IPATH"/*.count > "$countout"
  total=$(awk '{sum+=$2} END {print sum}' "$countout")
  echo "Wrote: $countout (total tests: $total)"
fi
