#!/bin/sh
set -e -u

blockMesh
rm -rfv 0/ && cp -r 0.orig/ 0/
touch solid-openfoam.foam

../../tools/run-openfoam.sh "$@"