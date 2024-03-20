#!/bin/bash
set -e -u

. ../../tools/log.sh

cp -r constant/polyMesh.orig constant/polyMesh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
