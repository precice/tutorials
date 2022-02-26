#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -fv ./*.vtk
clean_precice_logs .
