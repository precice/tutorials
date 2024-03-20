#!/bin/sh
set -e -u

. ../../tools/log.sh

log blockMesh

log ../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && log openfoam_remove_empty_dirs

close_log
