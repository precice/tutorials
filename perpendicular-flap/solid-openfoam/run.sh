#!/usr/bin/env bash

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

echo "[WARNING] This case is only provided for quick testing, it gives significantly different results than other cases in this tutorial, and does not work in parallel. See solid-solids4foam for a better OpenFOAM-based solid case.  Discussion on https://github.com/precice/tutorials/issues/515."
sleep 5

blockMesh

# Compile boundary condition
(cd solidDisplacementFoamForce && wmake libso)

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
