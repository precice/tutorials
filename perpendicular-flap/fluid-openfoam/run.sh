#!/bin/sh
set -e -u

blockMesh
touch fluid-openfoam.foam

../../tools/run-openfoam.sh "$@"
