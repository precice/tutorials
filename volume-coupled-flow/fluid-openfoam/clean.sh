#!/usr/bin/env sh
set -e -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

. ../../tools/cleaning-tools.sh

clean_openfoam .

# Necessary signal file for the system tests,
# to keep the Source container running till the end
rm -f ../fluid-participant-finished.log