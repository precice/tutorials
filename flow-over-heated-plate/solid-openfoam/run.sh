#!/bin/sh
set -e -u

blockMesh
touch solid-openfoam.foam

../../tools/run-openfoam.sh "$@"
