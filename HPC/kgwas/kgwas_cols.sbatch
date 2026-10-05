#!/bin/bash
#SBATCH --job-name=kgwas_cols
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G

set -euo pipefail

"${REPO_ROOT}/HPC/kgwas/kgwas_cols.sh"
