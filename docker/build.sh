#!/bin/bash
set -euo pipefail

IMAGE_NAME="${1:-onewaykgwas}"

# Change to repository root (parent of docker/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

SIF="${REPO_ROOT}/${IMAGE_NAME}.sif"

echo "[*] Building Docker image: ${IMAGE_NAME}:latest"
docker build -t "${IMAGE_NAME}:latest" -f docker/Dockerfile .

# Check for singularity/apptainer
if command -v apptainer >/dev/null 2>&1; then
  BUILDER="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  BUILDER="singularity"
else
  echo "[*] Done (Docker only - no apptainer/singularity found)"
  exit 0
fi

echo "[*] Building SIF with ${BUILDER}"
if ${BUILDER} build --force "${SIF}" "docker-daemon://${IMAGE_NAME}:latest" 2>/dev/null; then
  echo "[*] Built ${SIF} from docker-daemon://"
else
  echo "[*] docker-daemon:// unavailable; falling back to docker-archive://"
  TAR="${REPO_ROOT}/${IMAGE_NAME}.tar"
  docker save -o "${TAR}" "${IMAGE_NAME}:latest"
  ${BUILDER} build --force "${SIF}" "docker-archive://${TAR}"
  rm -f "${TAR}"
fi

echo "[*] Done: ${SIF}"

