#!/bin/sh
set -e -u

blockMesh
touch openfoam-neumann.foam
./setInitialField.sh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
