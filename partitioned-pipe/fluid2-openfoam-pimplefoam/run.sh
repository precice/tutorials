#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

use_skewed=false

for arg in "$@"
do
  shift
  [ "$arg" = "-skewed" ] && use_skewed=true && continue
  set -- "$@" "$arg"
done

if [ "$use_skewed" = true ]; then
    blockMesh -dict system/blockMeshDictSkewed
else
    blockMesh
fi

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

postProcess -func "flowRatePatch(name=inlet)" -latestTime -noZero

close_log
