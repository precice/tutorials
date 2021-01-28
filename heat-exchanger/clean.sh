#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

echo "Cleaning..."

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

# Participant 1: Inner-Fluid
Participant1="inner-fluid-openfoam"
cd ${Participant1}
    # Clean the case
    cleanCase
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
    rm -fv \
        precice-*.log \
        precice-postProcessingInfo.log \
        precice-*-events.json
cd ..

# Participant 2: Outer-Fluid
Participant2="outer-fluid-openfoam"
cd ${Participant2}
    # Clean the case
    cleanCase
    # Create an empty .foam file for ParaView
    # Note: ".foam" triggers the native OpenFOAM reader of ParaView.
    # Change to ".OpenFOAM" to use the OpenFOAM reader provided with OpenFOAM.
    touch ${Participant2}.foam
    # Remove the log files
    rm -fv ${Participant2}_blockMesh.log
    rm -fv ${Participant2}_checkMesh.log
    rm -fv ${Participant2}_decomposePar.log
    rm -fv ${Participant2}.log
    rm -fv ${Participant2}_reconstructPar.log
    rm -fv \
        precice-*.log \
        precice-postProcessingInfo.log \
        precice-*-events.json
cd ..

# Participant 3: Solid
Participant3="solid-calculix"
cd ${Participant3}
    # Delete result and log files
    rm -fv solid.cvg solid.dat solid.frd solid.sta
    rm -fv \
        precice-*.log \
        precice-postProcessingInfo.log \
        precice-*-events.json
cd ..

# Remove the preCICE address files
rm -rfv precice-run
rm -fv .*.address

echo "Cleaning complete!"
#------------------------------------------------------------------------------
