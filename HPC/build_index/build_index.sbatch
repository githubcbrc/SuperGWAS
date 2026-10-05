#!/bin/bash
#SBATCH --job-name=build_index
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=80G

set -euo pipefail

"${REPO_ROOT}/HPC/build_index/build_index.sh"
