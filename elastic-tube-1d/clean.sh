#!/bin/sh
set -e -u

echo "Cleaning logs and output files from previous runs..."

rm -f precice-Fluid-events.json \
      precice-Solid-events.json \
      precice-Fluid-iterations.log \
      precice-Solid-convergence.log \
      precice-Solid-iterations.log \
      precice-Fluid-events-summary.log \
      precice-STRUCTUR-events-summary.log \
      postproc/*.vtk \

rm -r precice-run/
rm -rfv precice-Solid-events-summary.log

echo "Cleaning successful!"
