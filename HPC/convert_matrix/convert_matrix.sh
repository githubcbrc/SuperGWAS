#!/usr/bin/env bash
#SBATCH --job-name=convert_matrix
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G

set -euo pipefail

# ---------------------------------------------------------
# Process ONE batch chosen by SLURM_ARRAY_TASK_ID.
# Direct binary invocation, no container.
# ---------------------------------------------------------

ts(){ date '+%Y-%m-%d %H:%M:%S'; }
log(){ echo "[$(ts)] $*"; }

if [[ -z "${REPO_ROOT:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi

if [[ -z "${ROW_MATRIX_DIR:-}" ]]; then
    source "${REPO_ROOT}/config.sh"
fi

BIN_DIR="${BIN_DIR:-${REPO_ROOT}/bin}"

# Settings
SORTED_OUT="${SORTED_OUT:-${RESULTS_ROOT}/sorted}"
BATCH_TMP="${BATCH_TMP:-${RESULTS_ROOT}/batches}"
MANIFEST="${SORTED_OUT}/MANIFEST"
JOBS_COMP="${JOBS_COMP:-$(nproc)}"
DRY_RUN="${DRY_RUN:-0}"

[[ -f "$MANIFEST" ]] || { echo "MANIFEST not found: $MANIFEST" >&2; exit 1; }

# Pick batch and validate
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
BNAME="$(sed -n "$((TASK_ID+1))p" "$MANIFEST" | tr -d '\r' | xargs -r)"
[[ -n "$BNAME" ]] || { echo "Invalid SLURM_ARRAY_TASK_ID=$TASK_ID for $MANIFEST" >&2; exit 1; }

BDIR="$BATCH_TMP/$BNAME"
OUTDIR="$SORTED_OUT/$BNAME"
BITSDIR="$OUTDIR/bits"

log "Task $TASK_ID -> $BNAME"
log "BDIR=$BDIR"
log "OUTDIR=$OUTDIR"

[[ -d "$BDIR" ]] || { echo "Missing batch dir: $BDIR" >&2; exit 1; }
mkdir -p "$OUTDIR"

# Resume: skip if the usual completion marker exists
if [[ -f "$OUTDIR/kmers.txt.gz" ]]; then
  log "Skip: $OUTDIR/kmers.txt.gz exists."
  exit 0
fi

# Compressor
if command -v pigz >/dev/null 2>&1; then COMP=(pigz -f -p "${JOBS_COMP}"); else COMP=(gzip -f); fi

if [[ "$DRY_RUN" == "1" ]]; then
  log "DRY_RUN=1 -> not executing."
  exit 0
fi

log "Executing convert_matrix..."
"${BIN_DIR}/convert_matrix" "$BDIR" "$OUTDIR"

# Post: compress outputs
if [[ -f "$OUTDIR/kmers.txt" ]]; then
  log "Compressing kmers.txt ..."
  "${COMP[@]}" "$OUTDIR/kmers.txt"
fi

if [[ -d "$BITSDIR" ]]; then
  log "Compressing bits/* ..."
  find "$BITSDIR" -maxdepth 1 -type f ! -name '*.gz' -print0 \
    | xargs -0 -r "${COMP[@]}"
fi

log "Done $BNAME"
