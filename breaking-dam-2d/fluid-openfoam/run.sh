#!/bin/sh
set -e -u

#cp -r 0_orig 0
#cp -r constant/polyMesh.orig constant/polyMesh
blockMesh
setFields
touch fluid-openfoam.foam

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
