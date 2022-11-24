#!/bin/bash

cd ${0%/*} || exit 1    		    		# Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions

echo "Preparing and running the Fluid2 participant..."

rm -rfv ./0/
cp -r ./0.orig/ ./0/
cd ..
blockMesh -case Fluid2
checkMesh -case Fluid2
cd Fluid2
# Run
solver=$(getApplication)

cd ..
$solver -case Fluid2 | tee Fluid2.log 2>&1

foamToVTK -case Fluid2

