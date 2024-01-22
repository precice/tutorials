#!/bin/sh
set -e -u

blockMesh
cp -r 0.orig 0
setFields
touch monolithic.foam

../../../tools/run-openfoam.sh "$@"
. ../../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
