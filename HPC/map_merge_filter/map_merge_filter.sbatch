#!/usr/bin/env bash
#SBATCH --job-name=map_merge_filter
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G

set -euo pipefail

"${REPO_ROOT}/HPC/map_merge_filter/map_merge_filter.sh"
