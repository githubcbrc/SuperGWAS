#!/usr/bin/env bash
# =============================================================================
# SuperGWAS 2W Two-Way Pipeline - submits all stages to SLURM with dependencies
# =============================================================================
#
# Usage:
#   ./run_pipeline_2w.sh [options]
#
# Options:
#   --interaction-only      Use 1-df interaction-only test (default: 3-df joint)
#   --skip-structure        Skip structure computation (use existing U0.txt)
#   --skip-projections      Skip projection computation (use existing Th/Tp)
#   --dry-run               Show what would be submitted without submitting
#
# Pipeline steps:
#   1. compute_structure_2w    - Dual PCA + Kronecker → U0.txt
#   2. compute_projections_2w  - Factored T projection → Th.txt, Tp.txt
#   3. kgwas_2w                - All-pairs host x pathogen sweep
#
# =============================================================================

set -euo pipefail

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export REPO_ROOT
source "${REPO_ROOT}/config.sh"

# Defaults
SKIP_STRUCTURE=false
SKIP_PROJECTIONS=false
DRY_RUN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --interaction-only)
            INTERACTION_ONLY=true
            shift
            ;;
        --skip-structure)
            SKIP_STRUCTURE=true
            shift
            ;;
        --skip-projections)
            SKIP_PROJECTIONS=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            head -21 "$0" | tail -17
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

HPC_DIR="${REPO_ROOT}/HPC"

# Create log directories
mkdir -p "${LOGS_DIR}"/{compute_structure_2w,compute_projections_2w,kgwas_2w}

# Auto-detect completed steps
if [[ "$SKIP_STRUCTURE" == "false" && -f "$U_2W" ]]; then
    echo "NOTE: Found existing ${U_2W}"
    echo "      Auto-skipping structure computation."
    echo ""
    SKIP_STRUCTURE=true
fi

if [[ "$SKIP_PROJECTIONS" == "false" && -f "$TH_FILE" && -f "$TP_FILE" ]]; then
    echo "NOTE: Found existing ${TH_FILE} and ${TP_FILE}"
    echo "      Auto-skipping projection computation."
    echo ""
    SKIP_PROJECTIONS=true
fi

echo "=============================================="
echo "SuperGWAS 2W Two-Way Pipeline"
echo "=============================================="
echo "Host matrix:       ${HOST_MATRIX}"
echo "Pathogen matrix:   ${PATHOGEN_MATRIX}"
echo "Phenotype:         ${PHENOTYPE_2W}"
echo "Rank host:         ${RANK_HOST}"
echo "Rank pathogen:     ${RANK_PATHOGEN}"
echo "Interaction only:  ${INTERACTION_ONLY}"
echo "Results root:      ${RESULTS_ROOT}"
echo ""
echo "Steps:"
echo "  1. compute_structure_2w:    $([ "$SKIP_STRUCTURE" = true ] && echo SKIP || echo RUN)"
echo "  2. compute_projections_2w:  $([ "$SKIP_PROJECTIONS" = true ] && echo SKIP || echo RUN)"
echo "  3. kgwas_2w:                RUN"
echo ""
echo "Dry run:           ${DRY_RUN}"
echo "=============================================="
echo ""

# Track job IDs for dependencies
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

# =============================================================================
# Step 1: compute_structure_2w
# =============================================================================
if [[ "$SKIP_STRUCTURE" == "false" ]]; then
    echo "Step 1: compute_structure_2w"
    LAST_JOB=$(submit_job "compute_structure_2w" \
        "${HPC_DIR}/compute_structure_2w/compute_structure_2w.sbatch" \
        "$LAST_JOB")
    echo "  Submitted job: $LAST_JOB"
    echo ""
else
    echo "Step 1: compute_structure_2w [SKIPPED]"
    echo ""
fi

# =============================================================================
# Step 2: compute_projections_2w
# =============================================================================
if [[ "$SKIP_PROJECTIONS" == "false" ]]; then
    echo "Step 2: compute_projections_2w"
    LAST_JOB=$(submit_job "compute_projections_2w" \
        "${HPC_DIR}/compute_projections_2w/compute_projections_2w.sbatch" \
        "$LAST_JOB")
    echo "  Submitted job: $LAST_JOB"
    echo ""
else
    echo "Step 2: compute_projections_2w [SKIPPED]"
    echo ""
fi

# =============================================================================
# Step 3: kgwas_2w
# =============================================================================
echo "Step 3: kgwas_2w"
LAST_JOB=$(submit_job "kgwas_2w" \
    "${HPC_DIR}/kgwas_2w/kgwas_2w.sbatch" \
    "$LAST_JOB")
echo "  Submitted job: $LAST_JOB"
echo ""

# =============================================================================
# Summary
# =============================================================================
echo "=============================================="
echo "2W Pipeline submitted!"
echo "=============================================="
echo ""
echo "Monitor with:"
echo "  watch 'squeue --me -o \"%.10i %.9P %.20j %.8u %.2t %.10M %.6D %R\"'"
echo ""
echo "Logs in: ${LOGS_DIR}/"
