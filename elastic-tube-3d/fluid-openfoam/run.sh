#!/bin/sh
set -e -u

cp -r constant/polyMesh.orig constant/polyMesh
touch fluid-openfoam.foam

../../tools/run-openfoam.sh "$@"
