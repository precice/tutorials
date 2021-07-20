#!/bin/sh
set -e -u

blockMesh
touch openfoam-neumann.foam

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
