#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -f ./results/Fluid1D_*
rm -f ./probes.txt
rm -f ./final_fields.txt
clean_precice_logs .
clean_case_logs .