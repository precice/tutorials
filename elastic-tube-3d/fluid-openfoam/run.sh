#!/bin/sh
set -e -u

. ../../tools/log.sh

log cp -r constant/polyMesh.orig constant/polyMesh

log ../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && log openfoam_remove_empty_dirs

close_log
