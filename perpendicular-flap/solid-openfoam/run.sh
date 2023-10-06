#!/bin/bash

blockMesh
touch solid-openfoam.foam

# Compile boundary condition
(cd solidDisplacementFoamForce && wmake libso)

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
