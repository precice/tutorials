#!/usr/bin/env bash
set -euo pipefail

# ---- Config (you can override via env vars) ----
REPO_HTTPS_DEFAULT="https://bitbucket.org/davidscn/mercurydpm.git"
REPO_URL="${REPO_URL:-$REPO_HTTPS_DEFAULT}"          # or git@bitbucket.org:davidscn/mercurydpm.git
BRANCH="${BRANCH:-channel-transport-tutorial}"
SRC_DIR="${SRC_DIR:-$PWD/mercurydpm}"
BUILD_DIR="${BUILD_DIR:-$SRC_DIR/build}"

# ---- Clone (if needed) and checkout branch ----
if [[ ! -d "$SRC_DIR/.git" ]]; then
  echo "Cloning repo into: $SRC_DIR"
  git clone "$REPO_URL" "$SRC_DIR"
else
  echo "Repo already exists at: $SRC_DIR"
fi

cd "$SRC_DIR"
git fetch --all --prune
git checkout "$BRANCH"
git pull --ff-only || true

# ---- Configure with CMake in build directory ----
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

echo "Configuring in: $BUILD_DIR"
cmake -D MercuryDPM_PreCICE_COUPLING="ON" ..

# ---- Build only the requested target ----
echo "Building target: ChannelTransport"
make -j "$(nproc)" ChannelTransport

export MERCURYDPM_BUILD_DIR="$BUILD_DIR"
echo "Exported MERCURYDPM_BUILD_DIR=$MERCURYDPM_BUILD_DIR"
echo "Done."

