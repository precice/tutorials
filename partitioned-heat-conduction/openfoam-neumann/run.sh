#!/bin/sh
set -e -u

blockMesh
topoSet
touch openfoam-neumann.foam

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
