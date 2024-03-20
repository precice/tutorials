#!/bin/sh
set -e -u

. ../../tools/log.sh

use_skewed=false

for arg in "$@"
do
  shift
  [ "$arg" = "-skewed" ] && use_skewed=true && continue
  set -- "$@" "$arg"
done

if [ "$use_skewed" = true ]; then
    log blockMesh -dict system/blockMeshDictSkewed
else
    log blockMesh
fi

log ../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && log openfoam_remove_empty_dirs

log postProcess -func "flowRatePatch(name=inlet)" -latestTime -noZero

close_log
