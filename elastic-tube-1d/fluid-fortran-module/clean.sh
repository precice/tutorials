#!/usr/bin/env sh
set -e -u

cd "$(cd "$(dirname "$0")" && pwd)"
. ../../tools/cleaning-tools.sh

rm -rvf ./output/*.vtk
clean_precice_logs .
clean_case_logs .
