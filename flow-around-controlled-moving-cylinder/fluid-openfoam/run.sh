#!/bin/sh
set -e -u

. ../../tools/log.sh

log blockMesh
log transformPoints -scale '(0.0016  0.0016 1)'
log transformPoints -translate '(0.0  0.0 -0.05)'

rm -rf 0
cp -r 0.orig 0

log ../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && log openfoam_remove_empty_dirs

close_log
