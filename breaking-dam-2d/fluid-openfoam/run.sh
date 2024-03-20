#!/bin/sh
set -e -u

cp 0/alpha.water_orig 0/alpha.water
blockMesh
setFields
touch fluid-openfoam.foam

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
