#!/usr/bin/env sh
set -e -u

# Determine the tutorial directory (where the symlink is located)
TUTORIAL_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$TUTORIAL_DIR"

# shellcheck disable=SC1091
. ../tools/cleaning-tools.sh

clean_tutorial .
clean_precice_logs .
rm -fv ./*.log
rm -fv ./*.vtu

