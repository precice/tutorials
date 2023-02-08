#!/bin/sh
set -e -u

if [ $# -eq 0 ]; then
    blockMesh 
elif [ "$1" = "skewed" ]; then
    blockMesh -dict system/blockMeshDictSkewed
else
    echo Error: Use \"./run.sh skewed\" to run the case with a skewed mesh
    exit 1
fi
touch fluid2-openfoam-pimplefoam.foam

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

postProcess -func "flowRatePatch(name=inlet)" -latestTime -noZero