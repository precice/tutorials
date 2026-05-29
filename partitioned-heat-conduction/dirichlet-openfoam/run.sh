#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

# Build the heatTransfer solver, if not found.
# Race condition: when running both participants with OpenFOAM,
# they both need the same solver to be built.
# The first one builds it, the second one waits.
if ! command -v heatTransfer > /dev/null 2>&1; then
  if [ ! -d "../solver-openfoam/Make/${WM_OPTIONS}" ]; then
    echo "Building the heatTransfer OpenFOAM solver"
    wmake ../solver-openfoam/
  else
    echo "The executable heatTransfer is not found, but the build directory ../solver-openfoam/Make/${WM_OPTIONS} was detected. A build is probably in progress, waiting 20s..."
    sleep 20
  fi
fi

blockMesh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
