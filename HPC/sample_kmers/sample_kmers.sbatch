#!/bin/bash
#SBATCH --job-name=sample_kmers
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

"${REPO_ROOT}/HPC/sample_kmers/sample_kmers.sh"
