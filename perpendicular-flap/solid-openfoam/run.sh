#!/bin/bash

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

blockMesh

# Compile boundary condition
(cd solidDisplacementFoamForce && wmake libso)

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
