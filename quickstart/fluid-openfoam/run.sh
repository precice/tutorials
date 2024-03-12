#!/bin/sh
set -e -u

. ../../tools/get-case-name.sh

blockMesh | tee "$CASENAME.log" 2>&1
touch "$CASENAME.foam"

../../tools/run-openfoam.sh "$@"  | tee --append "$CASENAME.log" 2>&1
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs | tee --append "$CASENAME.log" 2>&1
