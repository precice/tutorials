#!/bin/sh
set -e -u

blockMesh
transformPoints -scale '(0.0016  0.0016 1)'
transformPoints -translate '(0.0  0.0 -0.05)'

rm -rf 0
cp -r 0.orig 0

touch fluid-openfoam.foam

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
