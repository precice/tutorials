#!/bin/bash

. ../../tools/log.sh

log blockMesh

# Compile boundary condition
(cd solidDisplacementFoamForce && log wmake libso)

log ../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && log openfoam_remove_empty_dirs

close_log
