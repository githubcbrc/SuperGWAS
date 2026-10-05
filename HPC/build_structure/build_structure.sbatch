#!/bin/bash
#SBATCH --job-name=build_structure
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

"${REPO_ROOT}/HPC/build_structure/build_structure.sh"
