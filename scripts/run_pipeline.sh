#!/usr/bin/env bash
# =============================================================================
# SuperGWAS Full Pipeline (1D + 2D) - submits all stages to SLURM
# =============================================================================
#
# Usage:
#   ./run_pipeline.sh [options]
#
# Options:
#   --mode <cols|rows>      Use kgwas --mode cols or --mode rows (default: cols)
#   --shards-per-job N      Shards per sampling job (default: 50)
#   --skip-sampling         Skip sampling (use existing sampled kmers)
#   --no-2d                 Skip 2D epistasis steps (6-8)
#   --dry-run               Show what would be submitted without submitting
#
# Pipeline steps:
#   1. sample_kmers          - Reservoir sampling from row matrix
#   2. build_structure       - filter_kmers + compute_structure_all → U.bin + U2D.bin
#   3. kgwas                 - Run GWAS (produces 1D .tsv + 2D .2d.tsv)
#   --- 1D path ---
#   4. build_index           - Build 1D k-mer lookup index (compact, no bitstrings)
#   5. map_merge_filter      - Map to reference + density filter
#   --- 2D path (parallel with 4-5 after step 3) ---
#   6. build_index_2d        - Build small 2D index from .2d.tsv files
#   7. kgwas_2d              - All-pairs epistasis (array job) + merge
#   8. map_2d_pairs          - Translate k-mer pairs to genome positions per ref
#   --- Final ---
#   9. plot                  - 1D Manhattan + 2D 3D Manhattan
#
# =============================================================================

set -euo pipefail

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export REPO_ROOT
source "${REPO_ROOT}/config.sh"

# Defaults
MODE="cols"
SHARDS_PER_JOB="${SHARDS_PER_JOB:-50}"
SKIP_SAMPLING=false
NO_2D=false
DRY_RUN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)          MODE="$2"; shift 2 ;;
        --shards-per-job) SHARDS_PER_JOB="$2"; shift 2 ;;
        --skip-sampling) SKIP_SAMPLING=true; shift ;;
        --no-2d)         NO_2D=true; shift ;;
        --dry-run)       DRY_RUN=true; shift ;;
        -h|--help)       head -32 "$0" | tail -28; exit 0 ;;
        *)               echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ "$MODE" != "cols" && "$MODE" != "rows" ]]; then
    echo "ERROR: --mode must be 'cols' or 'rows'" >&2
    exit 1
fi

HPC_DIR="${REPO_ROOT}/HPC"

# Create log directories
LOG_DIRS=(sample_kmers build_structure kgwas build_index map_merge_filter plot)
if [[ "$NO_2D" == "false" ]]; then
    LOG_DIRS+=(build_index_2d kgwas_2d map_2d_pairs)
fi
for d in "${LOG_DIRS[@]}"; do mkdir -p "${LOGS_DIR}/${d}"; done

# =========================================================================
# Auto-detect skips
# =========================================================================
SKIP_STRUCTURE=false
SKIP_KGWAS=false
SKIP_BUILD_INDEX=false
SKIP_MAP_1D=false
SKIP_BUILD_INDEX_2D=false
SKIP_KGWAS_2D=false
SKIP_MAP_2D_PAIRS=false

if [[ "$SKIP_SAMPLING" == "false" && -f "$SAMPLED_KMERS" ]]; then
    echo "NOTE: Found existing ${SAMPLED_KMERS} — skipping sampling"
    SKIP_SAMPLING=true
fi

if [[ -f "$U_BIN" && -f "$U2D_BIN" ]]; then
    echo "NOTE: Found existing ${U_BIN} + ${U2D_BIN} — skipping build_structure"
    SKIP_STRUCTURE=true
fi

KGWAS_BATCH_DIR="${KGWAS_OUT}/${PHENO_NAME}"
if [[ -d "$KGWAS_BATCH_DIR" ]]; then
    KGWAS_TSV_COUNT=$(find "$KGWAS_BATCH_DIR" -maxdepth 1 -name '*.tsv' ! -name '*.2d.tsv' 2>/dev/null | wc -l)
    if [[ "$KGWAS_TSV_COUNT" -gt 0 ]]; then
        echo "NOTE: Found ${KGWAS_TSV_COUNT} .tsv files in ${KGWAS_BATCH_DIR} — skipping kgwas"
        SKIP_KGWAS=true
    fi
fi

INDEX_PATH="${INDEX_DIR}/${PHENO_NAME}"
if [[ -d "$INDEX_PATH" && -f "${INDEX_PATH}/index.meta" ]]; then
    echo "NOTE: Found existing ${INDEX_PATH}/index.meta — skipping build_index"
    SKIP_BUILD_INDEX=true
fi

MAPPING_1D_DIR="${MAPPING_OUT}/${PHENO_NAME}"
FILTERED_COUNT=$(find "$MAPPING_1D_DIR" -name '*.filtered.tsv' 2>/dev/null | wc -l || true)
if [[ "$FILTERED_COUNT" -gt 0 ]]; then
    echo "NOTE: Found ${FILTERED_COUNT} .filtered.tsv in ${MAPPING_1D_DIR} — skipping map_merge_filter"
    SKIP_MAP_1D=true
fi

INDEX_2D_DIR="${INDEX_DIR}/${PHENO_NAME}_2d"
export INDEX_2D_DIR
if [[ "$NO_2D" == "false" ]]; then
    if [[ -d "$INDEX_2D_DIR" && -f "${INDEX_2D_DIR}/index.meta" ]]; then
        echo "NOTE: Found existing ${INDEX_2D_DIR}/index.meta — skipping build_index_2d"
        SKIP_BUILD_INDEX_2D=true
    fi
    KGWAS_2D_TSV="${KGWAS_2D_OUT}/${PHENO_NAME}_2d.tsv"
    KGWAS_2D_KMERS="${KGWAS_2D_OUT}/${PHENO_NAME}_2d.kmers"
    if [[ -f "$KGWAS_2D_TSV" && -f "$KGWAS_2D_KMERS" ]]; then
        echo "NOTE: Found existing ${KGWAS_2D_TSV} + .kmers — skipping kgwas_2d"
        SKIP_KGWAS_2D=true
    fi
    SKIP_MAP_2D_PAIRS=false
    MAP_2D_COUNT=$(find "$KGWAS_2D_OUT" -mindepth 2 -name '*_2d.tsv' 2>/dev/null | wc -l || true)
    if [[ "$MAP_2D_COUNT" -gt 0 ]]; then
        echo "NOTE: Found ${MAP_2D_COUNT} per-ref mapped 2D pairs — skipping map_2d_pairs"
        SKIP_MAP_2D_PAIRS=true
    fi
fi

echo ""
echo "=============================================="
echo "SuperGWAS Pipeline (1D + 2D)"
echo "=============================================="
echo "Mode:            kgwas --mode ${MODE}"
echo "Phenotype:       ${PHENO_NAME}"
echo "Results root:    ${RESULTS_ROOT}"
echo ""
echo "Steps:"
echo "  1. sample_kmers:           $([ "$SKIP_SAMPLING" = true ] && echo SKIP || echo RUN)"
echo "  2. build_structure:        $([ "$SKIP_STRUCTURE" = true ] && echo SKIP || echo RUN)"
echo "  3. kgwas (${MODE}):            $([ "$SKIP_KGWAS" = true ] && echo SKIP || echo RUN)"
echo "  4. build_index:            $([ "$SKIP_BUILD_INDEX" = true ] && echo SKIP || echo RUN)"
echo "  5. map_merge_filter:       $([ "$SKIP_MAP_1D" = true ] && echo SKIP || echo RUN)"
if [[ "$NO_2D" == "false" ]]; then
    echo "  6. build_index_2d:         $([ "$SKIP_BUILD_INDEX_2D" = true ] && echo SKIP || echo RUN)"
    echo "  7. kgwas_2d (${NUM_2D_TASKS:-16} tasks): $([ "$SKIP_KGWAS_2D" = true ] && echo SKIP || echo RUN)"
    echo "  8. map_2d_pairs (per ref):  $([ "$SKIP_MAP_2D_PAIRS" = true ] && echo SKIP || echo RUN)"
else
    echo "  6-8. 2D steps:             DISABLED (--no-2d)"
fi
echo "  9. plot (per ref):         RUN"
echo ""
echo "Dry run:         ${DRY_RUN}"
echo "=============================================="
echo ""

# =========================================================================
# Job submission helper
# =========================================================================
LAST_JOB=""

submit_job() {
    local name="$1"
    local script="$2"
    local dep="${3:-}"
    local extra="${4:-}"

    local cmd="sbatch --parsable --export=ALL"
    cmd+=" -o ${LOGS_DIR}/${name}/%x.%j.out"
    cmd+=" -e ${LOGS_DIR}/${name}/%x.%j.err"

    if [[ -n "$dep" ]]; then
        cmd+=" --dependency=afterok:${dep}"
    fi

    if [[ -n "$extra" ]]; then
        cmd+=" ${extra}"
    fi

    cmd+=" ${script}"

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] $cmd" >&2
        echo "fake_job_${name}"
    else
        echo "[SUBMIT] $name" >&2
        eval "$cmd"
    fi
}

# =========================================================================
# Step 1: sample_kmers
# =========================================================================
if [[ "$SKIP_SAMPLING" == "false" ]]; then
    echo "Step 1: sample_kmers"
    NUM_SHARDS=$(ls -1 "${ROW_MATRIX_DIR}"/*.tsv.gz 2>/dev/null | wc -l)
    if [[ "$NUM_SHARDS" -eq 0 ]]; then
        if [[ "$DRY_RUN" == "true" ]]; then NUM_SHARDS=100
        else echo "ERROR: No .tsv.gz files in ${ROW_MATRIX_DIR}" >&2; exit 1; fi
    fi
    NUM_JOBS=$(( (NUM_SHARDS + SHARDS_PER_JOB - 1) / SHARDS_PER_JOB ))
    MAX_IDX=$((NUM_JOBS - 1))
    echo "  Shards: $NUM_SHARDS, Jobs: $NUM_JOBS"
    LAST_JOB=$(submit_job "sample_kmers" \
        "${HPC_DIR}/sample_kmers/sample_kmers.sbatch" \
        "$LAST_JOB" \
        "--array=0-${MAX_IDX} --export=ALL,SHARDS_PER_JOB=${SHARDS_PER_JOB}")
    echo "  Submitted: $LAST_JOB"
    echo ""
else
    echo "Step 1: sample_kmers [SKIPPED]"
    echo ""
fi

# =========================================================================
# Step 2: build_structure (filter_kmers + compute_structure_all → U.bin + U2D.bin)
# =========================================================================
if [[ "$SKIP_STRUCTURE" == "false" ]]; then
    echo "Step 2: build_structure"
    LAST_JOB=$(submit_job "build_structure" \
        "${HPC_DIR}/build_structure/build_structure.sbatch" \
        "$LAST_JOB")
    echo "  Submitted: $LAST_JOB"
    echo ""
else
    echo "Step 2: build_structure [SKIPPED]"
    echo ""
fi

# =========================================================================
# Step 3: kgwas (produces .tsv + .2d.tsv)
# =========================================================================
if [[ "$SKIP_KGWAS" == "false" ]]; then
    echo "Step 3: kgwas (${MODE})"
    if [[ "$MODE" == "cols" ]]; then
        NUM_BATCHES=$(set +o pipefail; ls -1d "${COL_MATRIX_DIR}"/batch_* 2>/dev/null | wc -l)
        if [[ "$NUM_BATCHES" -eq 0 ]]; then
            if [[ "$DRY_RUN" == "true" ]]; then NUM_BATCHES=10
            else echo "ERROR: No batches in ${COL_MATRIX_DIR}" >&2; exit 1; fi
        fi
        MAX_IDX=$((NUM_BATCHES - 1))
        KGWAS_JOB=$(submit_job "kgwas" \
            "${HPC_DIR}/kgwas/kgwas_cols.sbatch" \
            "$LAST_JOB" \
            "--array=0-${MAX_IDX}")
    else
        NUM_CHUNKS=$(set +o pipefail; ls -1 "${ROW_MATRIX_DIR}"/*.tsv* 2>/dev/null | wc -l)
        if [[ "$NUM_CHUNKS" -eq 0 ]]; then
            if [[ "$DRY_RUN" == "true" ]]; then NUM_CHUNKS=10
            else echo "ERROR: No chunks in ${ROW_MATRIX_DIR}" >&2; exit 1; fi
        fi
        MAX_IDX=$((NUM_CHUNKS - 1))
        KGWAS_JOB=$(submit_job "kgwas" \
            "${HPC_DIR}/kgwas/kgwas_rows.sbatch" \
            "$LAST_JOB" \
            "--array=0-${MAX_IDX}")
    fi
    echo "  Submitted: $KGWAS_JOB (0-${MAX_IDX})"
    LAST_JOB="$KGWAS_JOB"
    echo ""
else
    echo "Step 3: kgwas [SKIPPED]"
    echo ""
fi

# Save kgwas job — 1D and 2D paths both branch from here
KGWAS_JOB_ID="$LAST_JOB"

# =========================================================================
# 1D path: steps 4-5
# =========================================================================

# --- Step 4: build_index (1D) ---
if [[ "$SKIP_BUILD_INDEX" == "false" ]]; then
    echo "Step 4: build_index"
    LAST_JOB=$(submit_job "build_index" \
        "${HPC_DIR}/build_index/build_index.sbatch" \
        "$LAST_JOB")
    echo "  Submitted: $LAST_JOB"
    echo ""
else
    echo "Step 4: build_index [SKIPPED]"
    echo ""
fi

# --- Step 5: map_merge_filter ---
NUM_REFS=$(find "${REF_DIR}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | wc -l)
if [[ "$NUM_REFS" -eq 0 ]]; then
    if [[ "$DRY_RUN" == "true" ]]; then NUM_REFS=1
    else echo "ERROR: No .fasta/.fa files in ${REF_DIR}" >&2; exit 1; fi
fi
MAX_IDX=$((NUM_REFS - 1))

if [[ "$SKIP_MAP_1D" == "false" ]]; then
    echo "Step 5: map_merge_filter"
    echo "  References: $NUM_REFS"
    LAST_JOB=$(submit_job "map_merge_filter" \
        "${HPC_DIR}/map_merge_filter/map_merge_filter.sbatch" \
        "$LAST_JOB" \
        "--array=0-${MAX_IDX}")
    echo "  Submitted: $LAST_JOB (0-${MAX_IDX})"
    echo ""
else
    echo "Step 5: map_merge_filter [SKIPPED]"
    echo ""
fi

LAST_JOB_1D="$LAST_JOB"

# =========================================================================
# 2D path: steps 6-8 (parallel with 4-5, branches from kgwas)
# =========================================================================
LAST_JOB_2D=""

if [[ "$NO_2D" == "false" ]]; then

    # --- Step 6: build_index_2d ---
    LAST_JOB_2D="$KGWAS_JOB_ID"
    if [[ "$SKIP_BUILD_INDEX_2D" == "false" ]]; then
        echo "Step 6: build_index_2d"
        LAST_JOB_2D=$(submit_job "build_index_2d" \
            "${HPC_DIR}/build_index_2d/build_index_2d.sbatch" \
            "$LAST_JOB_2D")
        echo "  Submitted: $LAST_JOB_2D"
        echo ""
    else
        echo "Step 6: build_index_2d [SKIPPED]"
        echo ""
    fi

    # --- Step 7: kgwas_2d (array job in k-mer space + merge) ---
    NUM_2D_TASKS="${NUM_2D_TASKS:-16}"
    export NUM_2D_TASKS
    if [[ "$SKIP_KGWAS_2D" == "false" ]]; then
        echo "Step 7: kgwas_2d (${NUM_2D_TASKS} tasks)"
        MAX_2D_IDX=$((NUM_2D_TASKS - 1))
        LAST_JOB_2D=$(submit_job "kgwas_2d" \
            "${HPC_DIR}/kgwas_2d/kgwas_2d.sbatch" \
            "$LAST_JOB_2D" \
            "--array=0-${MAX_2D_IDX}")
        echo "  Submitted: $LAST_JOB_2D (0-${MAX_2D_IDX})"

        echo "Step 7 (merge):"
        LAST_JOB_2D=$(submit_job "kgwas_2d" \
            "${HPC_DIR}/kgwas_2d/kgwas_2d_merge.sbatch" \
            "$LAST_JOB_2D")
        echo "  Submitted: $LAST_JOB_2D"
        echo ""
    else
        echo "Step 7: kgwas_2d [SKIPPED]"
        echo ""
    fi

    # --- Step 8: map_2d_pairs (per ref, lightweight) ---
    if [[ "$SKIP_MAP_2D_PAIRS" == "false" ]]; then
        echo "Step 8: map_2d_pairs (per ref)"
        LAST_JOB_2D=$(submit_job "map_2d_pairs" \
            "${HPC_DIR}/map_2d_pairs/map_2d_pairs.sbatch" \
            "$LAST_JOB_2D" \
            "--array=0-${MAX_IDX}")
        echo "  Submitted: $LAST_JOB_2D (0-${MAX_IDX})"
        echo ""
    else
        echo "Step 8: map_2d_pairs [SKIPPED]"
        echo ""
    fi
fi

# =========================================================================
# Step 9: plot (1D + 2D per ref, depends on both 1D and 2D completion)
# =========================================================================
echo "Step 9: plot (per ref)"

# Build dependency: wait for both 1D and 2D paths
PLOT_DEP=""
[[ -n "$LAST_JOB_1D" ]] && PLOT_DEP="$LAST_JOB_1D"
if [[ "$NO_2D" == "false" && -n "$LAST_JOB_2D" && "$LAST_JOB_2D" != "$KGWAS_JOB_ID" ]]; then
    if [[ -n "$PLOT_DEP" ]]; then
        PLOT_DEP="${PLOT_DEP}:${LAST_JOB_2D}"
    else
        PLOT_DEP="$LAST_JOB_2D"
    fi
fi

PLOT_JOB=$(submit_job "plot" \
    "${HPC_DIR}/plot/plot.sbatch" \
    "$PLOT_DEP" \
    "--array=0-${MAX_IDX}")
echo "  Submitted: $PLOT_JOB (0-${MAX_IDX})"
echo ""

# =========================================================================
# Summary
# =========================================================================
echo "=============================================="
echo "Pipeline submitted!"
echo "=============================================="
echo ""
echo "Monitor with:"
echo "  watch 'squeue --me -o \"%.10i %.9P %.20j %.8u %.2t %.10M %.6D %R\"'"
echo ""
echo "Logs in: ${LOGS_DIR}/"
