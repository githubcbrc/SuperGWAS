#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

"${REPO_ROOT}/HPC/plot/plot.sh"
