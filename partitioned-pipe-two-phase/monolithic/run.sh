#!/bin/sh
set -e -u

. ../../tools/log.sh

log blockMesh
log cp -r 0.orig 0
log setFields

log ../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && log openfoam_remove_empty_dirs

close_log
