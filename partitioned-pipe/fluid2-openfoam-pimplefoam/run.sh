#!/bin/sh
set -e -u

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

touch fluid2-openfoam-pimplefoam.foam

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

postProcess -func "flowRatePatch(name=inlet)" -latestTime -noZero