#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

rm ./results/Fluid1D_*
clean_precice_logs .
clean_case_logs .