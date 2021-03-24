#!/bin/bash
cd ${0%/*} || exit 1    		    # Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions

# Fluid participant

# Run this script in one terminal and the execute the runSolid1 and runSolid2 in separate terminals.
# These scripts present how the three participants would be started manually.

# Run this script with "-parallel" for parallel simulations

# The script "Allclean" cleans-up the result and log files.

# 1 for true, 0 for false
parallel=0
if [ "$1" = "-parallel" ]; then
    parallel=1
fi

echo "Preparing and running the Fluid participant..."

rm -rfv 0/
cp -r 0.orig/ 0/
blockMesh 
checkMesh

# Run
cd Fluid
	solver=$(getApplication)
	procs=$(getNumberOfProcessors)
cd ..
if [ $parallel -eq 1 ]; then
    decomposePar -force
    mpirun -np $procs $solver -parallel
    reconstructPar
else
    $solver
fi

# Workaround for issue #26 (OF-adapter, relevant for OF .com versions) 
./removeObsoleteFolders.sh

echo ""
echo "### NOTE ### Make sure to use the correct solver for your OpenFOAM version! (pimpleFoam for OpenFOAM v1806, OpenFOAM 6, or newer, vs pimpleDyMFoam for older) You may change this in your Fluid/system/controlDict file, if needed."
