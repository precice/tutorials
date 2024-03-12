#!/bin/sh
set -e -u

. ../../tools/log.sh

log blockMesh
touch "$CASENAME.foam"

log ../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

log date
