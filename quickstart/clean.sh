#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

echo "Cleaning..."

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

# Participant 1: Fluid (OpenFOAM)
Participant1="Fluid"
cd fluid-openfoam
    # Clean the case
    echo "Cleaning the fluid-openfoam case"
    cleanCase
    rm -rfv 0
    # Create an empty .foam file for ParaView
    # Note: ".foam" triggers the native OpenFOAM reader of ParaView.
    # Change to ".OpenFOAM" to use the OpenFOAM reader provided with OpenFOAM.
    touch ${Participant1}.foam

    # Remove the log files
    rm -fv ${Participant1}_blockMesh.log
    rm -fv ${Participant1}_checkMesh.log
    rm -fv ${Participant1}_decomposePar.log
    rm -fv ${Participant1}.log
    rm -fv ${Participant1}_reconstructPar.log
    
    # Remove the preCICE-related log files
    rm -fv \
    precice-*.log \
    precice-*-events.json
    
cd ..

# Participant 2: Solid (rigid body)
Participant2="Solid"
cd  solid-cpp
    # Clean the case
    echo "Cleaning the solid-cpp case"
    
    # Remove the preCICE-related log files
    echo "Deleting the preCICE log files..."
    rm -fv \
    precice-*.log \
    precice-*-events.json

    # Output files for preCICE versions before 1.2:
    rm -fv \
    iterations-${Participant1}.txt iterations-${Participant2}.txt \
    convergence-${Participant1}.txt convergence-${Participant2}.txt \
    Events-${Participant1}.log Events-${Participant2}.log \
    EventTimings-${Participant1}.log EventTimings-${Participant2}.log

cd ..

rm -rfv precice-run
rm -fv .${Participant1}-${Participant2}.address
rm -rf coupling-meshes

echo "Cleaning complete!"
#------------------------------------------------------------------------------
