#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -rvf ./output/*.vtk
clean_precice_logs .
