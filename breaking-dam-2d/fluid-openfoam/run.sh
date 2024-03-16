#!/bin/sh
set -e -u

. ../../tools/log.sh

log cp 0/alpha.water_orig 0/alpha.water
log blockMesh
log setFields

log ../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && log openfoam_remove_empty_dirs

close_log
