#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -f ./results/Fluid1D_*
rm -f ./probes.txt
clean_precice_logs .
clean_case_logs .
