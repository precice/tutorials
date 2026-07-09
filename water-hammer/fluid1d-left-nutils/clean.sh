#!/usr/bin/env sh
set -e -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

. ../../tools/cleaning-tools.sh

rm -f ./results/Fluid1D_*
rm -f ./probes.txt
clean_precice_logs .
clean_case_logs .
