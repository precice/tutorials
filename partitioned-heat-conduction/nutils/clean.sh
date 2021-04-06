#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -f Dirichlet-*.vtk
rm -f Neumann-*.vtk
rm -f *.log
rm -f *-events.json
