#!/bin/bash
#SBATCH --job-name=plot_1d
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

"${REPO_ROOT}/HPC/plot/plot_1d.sh"
