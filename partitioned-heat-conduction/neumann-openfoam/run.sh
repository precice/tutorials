#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if ! command -v heatTransfer > /dev/null 2>&1; then
  echo "Building the heatTransfer OpenFOAM solver"
  wclean ../solver-openfoam/
  wmake ../solver-openfoam/
fi

blockMesh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
