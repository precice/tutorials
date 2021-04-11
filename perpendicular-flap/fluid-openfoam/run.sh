#!/bin/sh
set -e -u

blockMesh
touch fluid.foam

../../tools/run-openfoam.sh "$@"
