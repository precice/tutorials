#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

blockMesh
transformPoints -scale '(0.0016  0.0016 1)'
transformPoints -translate '(0.0  0.0 -0.05)'

rm -rf 0
cp -r 0.orig 0

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
