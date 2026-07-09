#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

# Build the dynamicScalarTransportFoam solver, if not found.
if ! command -v dynamicScalarTransportFoam > /dev/null 2>&1; then
  echo "Building ../solver-openfoam and installing to ${FOAM_USER_APPBIN}/dynamicScalarTransportFoam..."
  wmake ../solver-openfoam
fi

blockMesh
cp -r 0.orig 0
setExprFields -time 0

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
