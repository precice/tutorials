#!/usr/bin/env sh
set -e -u

# Always operate relative to the directory of this script
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "Cleaning up all tutorials..."

find . -maxdepth 2 -mindepth 2 -name clean-tutorial.sh -execdir sh -c './clean-tutorial.sh' \;
