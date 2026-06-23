#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

# Build the heatTransfer solver, if not found.
if ! command -v heatTransfer > /dev/null 2>&1; then
  echo "Building ../solver-openfoam in a temporary build directory and installing to ${FOAM_USER_APPBIN}/heatTransfer..."
  cp -r ../solver-openfoam ./_solver-openfoam-copy
  wmake ./_solver-openfoam-copy
  rm -r ./_solver-openfoam-copy
fi

blockMesh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
