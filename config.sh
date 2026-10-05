#!/usr/bin/env bash
# =============================================================================
# SuperGWAS Configuration File
# =============================================================================
# Copy this file to config.local.sh and edit paths for your environment.
# Scripts will source config.local.sh if it exists, otherwise config.sh.
#
# Input paths are independent - set each to wherever your data lives.
# Output paths are derived from RESULTS_ROOT.
# =============================================================================

# If config.local.sh exists, source it instead and return
_CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${_CONFIG_DIR}/config.local.sh" && "${BASH_SOURCE[0]}" != *"config.local.sh" ]]; then
    source "${_CONFIG_DIR}/config.local.sh"
    return 0 2>/dev/null || exit 0
fi

# -----------------------------------------------------------------------------
# Container Settings
# -----------------------------------------------------------------------------
RUNNER="singularity"                    # "singularity" or "docker"
IMAGE_SIF="/path/to/onewaykgwas.sif"    # Singularity image path
IMAGE_DOCKER="onewaykgwas:latest"       # Docker image name

# -----------------------------------------------------------------------------
# Input Files (absolute paths - can be anywhere on filesystem)
# -----------------------------------------------------------------------------
ACCESSIONS="/path/to/samples/all_accessions.txt"
PHENOTYPE="/path/to/phenotypes/trichome_phenotype.txt"
REF_DIR="/path/to/refs"                         # Directory containing reference .fasta files

# -----------------------------------------------------------------------------
# Input Directories (absolute paths - can be anywhere on filesystem)
# -----------------------------------------------------------------------------
ROW_MATRIX_DIR="/path/to/kmer_matrix_rows"  # Row-format matrix (for sampling and kgwas_rows)
COL_MATRIX_DIR="/path/to/kmer_matrix_cols"  # Column-format matrix (for kgwas_cols)

# -----------------------------------------------------------------------------
# Output Root (all outputs go under this directory)
# -----------------------------------------------------------------------------
RESULTS_ROOT="/path/to/results"

# Derived output directories and files
STRUCTURE_OUT="${RESULTS_ROOT}/structure"
U_BIN="${STRUCTURE_OUT}/U.bin"
SAMPLED_OUT="${RESULTS_ROOT}/sampled"
SAMPLED_KMERS="${SAMPLED_OUT}/sampled_kmers.tsv.gz"
KGWAS_OUT="${RESULTS_ROOT}/kgwas"
MAPPING_OUT="${RESULTS_ROOT}/mappings"
INDEX_DIR="${RESULTS_ROOT}/index"
LOGS_DIR="${RESULTS_ROOT}/logs"

# -----------------------------------------------------------------------------
# Analysis Parameters
# -----------------------------------------------------------------------------
PHENO_NAME="trichome"                   # Phenotype name (used in output files)
K=51                                    # K-mer size
RANK=3                                  # PCA rank for compute_structure
LOG10_THRESH=6                          # -log10(p) threshold for output
CORR_THRESHOLD=0.2                      # Correlation pre-filter (0 = disabled)

# -----------------------------------------------------------------------------
# Mapping Parameters
# -----------------------------------------------------------------------------
LOG10_MIN=6                             # Min -log10(p) for mapping

# -----------------------------------------------------------------------------
# Sampling Parameters
# -----------------------------------------------------------------------------
RESERVOIR_SIZE=10000000                 # K-mer reservoir sample size
SUBSAMPLE_SIZE=200000                   # Matrix subsample size
MAF_THRESHOLD=5                         # Minor allele frequency threshold
SEED=42                                 # Random seed
SHARDS_PER_JOB=50                       # Shards per job (distributed sampling)

# -----------------------------------------------------------------------------
# Threading
# -----------------------------------------------------------------------------
THREADS=48                              # General thread count
LOAD_THREADS=24                         # Threads for data loading
BLAS_THREADS=112                        # Threads for BLAS operations
EMIT_THREADS=112                        # Threads for output writing

# -----------------------------------------------------------------------------
# Post-processing
# -----------------------------------------------------------------------------
WINDOW_SIZE=100                         # Density filter window size
MIN_DENSITY=0.95                        # Density filter threshold
CHUNK_SIZE=10000000                     # Intra-contig chunk size for parallelism (10 Mbp)

# -----------------------------------------------------------------------------
# 2D Epistasis Parameters
# -----------------------------------------------------------------------------
RANK_2D=3                               # PCA rank for 2D structure
LOG10_THRESH_2D_INPUT=15                # -log10(p) to select k-mers for 2D analysis
LOG10_THRESH_2D=18                      # -log10(p) threshold for kgwas_2d output
STRUCTURE_2D_OUT="${RESULTS_ROOT}/structure_2d"
U2D_BIN="${STRUCTURE_2D_OUT}/U2D.bin"
KGWAS_2D_OUT="${RESULTS_ROOT}/kgwas_2d"

# -----------------------------------------------------------------------------
# 2W Two-Way (host × pathogen) Parameters
# -----------------------------------------------------------------------------
HOST_MATRIX="/path/to/host_matrix"              # Host k-mer bitstring matrix
PATHOGEN_MATRIX="/path/to/pathogen_matrix"      # Pathogen k-mer bitstring matrix
PHENOTYPE_2W="/path/to/phenotype_2w.txt"        # Two-way phenotype file
RANK_HOST=3                             # PCA rank for host
RANK_PATHOGEN=3                         # PCA rank for pathogen
INTERACTION_ONLY=false                  # Use interaction-only test (1df vs 3df)
STRUCTURE_2W_OUT="${RESULTS_ROOT}/structure_2w"
U_2W="${STRUCTURE_2W_OUT}/U0.txt"
TH_FILE="${STRUCTURE_2W_OUT}/Th.txt"
TP_FILE="${STRUCTURE_2W_OUT}/Tp.txt"
KGWAS_2W_OUT="${RESULTS_ROOT}/kgwas_2w"

# =============================================================================
# Helper: Collect unique mount points from all paths
# =============================================================================
get_bind_paths() {
    local all_paths=(
        "$ACCESSIONS" "$PHENOTYPE" "$REF_DIR"
        "$ROW_MATRIX_DIR" "$COL_MATRIX_DIR" "$RESULTS_ROOT"
        "$HOST_MATRIX" "$PATHOGEN_MATRIX" "$PHENOTYPE_2W"
    )

    local -A mount_dirs
    for p in "${all_paths[@]}"; do
        if [[ -n "$p" && "$p" != "/path/to"* ]]; then
            # Use first 4 components as mount point (e.g., /scratch/user/data)
            # Walk up to find the deepest existing ancestor
            local mount_point="$p"
            while [[ -n "$mount_point" && ! -d "$mount_point" ]]; do
                mount_point="$(dirname "$mount_point")"
            done
            [[ -n "$mount_point" && "$mount_point" != "/" ]] && mount_dirs["$mount_point"]=1
        fi
    done

    local bind_str=""
    for dir in "${!mount_dirs[@]}"; do
        [[ -n "$bind_str" ]] && bind_str="${bind_str},"
        bind_str="${bind_str}${dir}:${dir}"
    done
    echo "$bind_str"
}

# =============================================================================
# Helper: Run command in container with auto-mounting
# =============================================================================
run_container() {
    local cmd="$1"
    shift
    local bind_str
    bind_str=$(get_bind_paths)

    local -a env_flags=(
        --env "CORR_THRESHOLD=${CORR_THRESHOLD}"
        --env "OMP_NUM_THREADS=${THREADS}"
        --env "OMP_PROC_BIND=spread"
        --env "OMP_PLACES=cores"
        --env "LOAD_THREADS=${LOAD_THREADS}"
        --env "BLAS_THREADS=${BLAS_THREADS}"
        --env "EMIT_THREADS=${EMIT_THREADS}"
    )

    if [[ "$RUNNER" == "singularity" ]]; then
        singularity exec --cleanenv \
            ${bind_str:+--bind "$bind_str"} \
            "${env_flags[@]}" \
            "$IMAGE_SIF" \
            "$cmd" "$@"
    else
        local -a docker_vols=()
        local -a MOUNTS
        IFS=',' read -ra MOUNTS <<< "$bind_str"
        for m in "${MOUNTS[@]}"; do
            docker_vols+=(-v "$m")
        done
        local -a docker_env=()
        for e in "${env_flags[@]}"; do
            # Convert --env "K=V" to -e "K=V"
            if [[ "$e" == "--env" ]]; then continue; fi
            docker_env+=(-e "$e")
        done
        docker run --rm "${docker_vols[@]}" \
            "${docker_env[@]}" \
            "$IMAGE_DOCKER" \
            "$cmd" "$@"
    fi
}

# Export all variables
export RUNNER IMAGE_SIF IMAGE_DOCKER
export ACCESSIONS PHENOTYPE REF_DIR
export ROW_MATRIX_DIR COL_MATRIX_DIR
export RESULTS_ROOT STRUCTURE_OUT U_BIN SAMPLED_OUT SAMPLED_KMERS KGWAS_OUT MAPPING_OUT INDEX_DIR LOGS_DIR
export PHENO_NAME K RANK LOG10_THRESH CORR_THRESHOLD
export LOG10_MIN
export RESERVOIR_SIZE SUBSAMPLE_SIZE MAF_THRESHOLD SEED SHARDS_PER_JOB
export THREADS LOAD_THREADS BLAS_THREADS EMIT_THREADS
export WINDOW_SIZE MIN_DENSITY CHUNK_SIZE
export RANK_2D LOG10_THRESH_2D_INPUT LOG10_THRESH_2D STRUCTURE_2D_OUT U2D_BIN KGWAS_2D_OUT
export HOST_MATRIX PATHOGEN_MATRIX PHENOTYPE_2W RANK_HOST RANK_PATHOGEN INTERACTION_ONLY
export STRUCTURE_2W_OUT U_2W TH_FILE TP_FILE KGWAS_2W_OUT
