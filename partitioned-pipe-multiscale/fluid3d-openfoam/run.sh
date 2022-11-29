#!/bin/bash

cd ${0%/*} || exit 1    		    		# Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions

echo "Preparing and running the Fluid3D participant..."

rm -rfv ./0/
cp -r ./0.orig/ ./0/
cd ..
blockMesh -case fluid3d-openfoam
checkMesh -case fluid3d-openfoam
cd fluid3d-openfoam
# Run
solver=$(getApplication)

cd ..
$solver -case fluid3d-openfoam | tee Fluid3D.log 2>&1

foamToVTK -case fluid3d-openfoam

