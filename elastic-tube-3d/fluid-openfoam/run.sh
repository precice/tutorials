#!/bin/sh
set -e -u

cp -r constant/polyMesh.orig constant/polyMesh
rm -rfv 0/ && cp -r 0.orig/ 0/
touch fluid-openfoam.foam

../../tools/run-openfoam.sh "$@"
