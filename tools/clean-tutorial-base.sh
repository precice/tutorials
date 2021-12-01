#!/bin/sh
set -e -u

# shellcheck disable=SC1091
. ../tools/cleaning-tools.sh

rm -fv result.vtk 
rm -fvr fine_mesh coarse_mesh mapped
clean_tutorial .
clean_precice_logs .
rm -fv *.log

