#!/bin/bash
set -e -u

. ../../tools/log.sh


blockMesh
./setInitialField.sh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
