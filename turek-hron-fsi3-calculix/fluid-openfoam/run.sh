#!/bin/sh
set -e -u

rm -rfv 0/
cp -r 0.orig/ 0/

blockMesh
checkMesh > Fluid_checkMesh.log
potentialFoam > Fluid_potentialFoam.log
touch fluid-openfoam.foam

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
