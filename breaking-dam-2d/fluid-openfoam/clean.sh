#!/usr/bin/env sh
set -e -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

. ../../tools/cleaning-tools.sh

clean_openfoam .
rm -fv 0/alpha.water
