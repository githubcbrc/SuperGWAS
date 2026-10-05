#!/bin/bash
#SBATCH --job-name=build_index_2d
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

set -euo pipefail

"${REPO_ROOT}/HPC/build_index_2d/build_index_2d.sh"
