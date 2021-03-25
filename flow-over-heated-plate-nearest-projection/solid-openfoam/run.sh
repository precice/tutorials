#!/bin/bash
cd ${0%/*} || exit 1    		    		# Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions

# Fluid participant

# Run this script in one terminal in order to start this participant.
# Run this script with "-parallel" for parallel simulations

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
touch solid-openfoam.foam

# Run
solver=$(getApplication)
procs=$(getNumberOfProcessors)

if [ $parallel -eq 1 ]; then
    decomposePar -force
    mpirun -np $procs $solver -parallel
    reconstructPar
else
    $solver
fi

# Workaround for issue #26
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

echo ""
echo "### NOTE ### Make sure to use the correct solver for your OpenFOAM version! (pimpleFoam for OpenFOAM v1806, OpenFOAM 6, or newer, vs pimpleDyMFoam for older) You may change this in your Fluid/system/controlDict file, if needed."
