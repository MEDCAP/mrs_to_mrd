#!/usr/bin/env bash
set -euo pipefail

# Ensure we have git commit hash
COMMIT_HASH=$(git rev-parse --short HEAD)
REPO="ghcr.io/medcap"
TARGETS=("convert" "shift" "recon")

# Check for uncommitted changes to avoid pushing stale code without warning
if ! git diff --quiet; then
  echo "WARNING: Working tree has uncommitted modifications."
fi

echo "================================================="
echo "Building and pushing for commit: ${COMMIT_HASH}"
echo "================================================="

for target in "${TARGETS[@]}"; do
  image="${REPO}/mrs-${target}"
  tag_commit="${image}:${COMMIT_HASH}"
  tag_latest="${image}:latest"

  echo ""
  echo "===== Building mrs-${target} (linux/amd64) ====="
  docker build \
    --platform linux/amd64 \
    --target "${target}" \
    -t "${tag_commit}" \
    -t "${tag_latest}" \
    .

  echo "===== Pushing mrs-${target} ====="
  docker push "${tag_commit}" 2>&1 | tail -n 4
  docker push "${tag_latest}" 2>&1 | tail -n 4
done

echo ""
echo "Done: All images successfully built and pushed."
