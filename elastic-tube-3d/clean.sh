#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

echo "Cleaning..."

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

# Participant 1: Fluid
Participant1="fluid-openfoam"
cd ${Participant1}
    # Clean the case
    # prevent cleaning the mesh
    mv ./constant/polyMesh ./constant/mesh
    cleanCase
    rm -rfv 0
    # restore mesh files
    mv ./constant/mesh ./constant/polyMesh
    # Create an empty .foam file for ParaView
    # Note: ".foam" triggers the native OpenFOAM reader of ParaView.
    # Change to ".OpenFOAM" to use the OpenFOAM reader provided with OpenFOAM.
    touch ${Participant1}.foam
    rm -fv ${Participant1}_decomposePar.log
    rm -fv ${Participant1}.log
    rm -fv ${Participant1}_reconstructPar.log
    rm -fv precice-*.log
    rm -fv precice-*-events.json
cd ..

# Participant 2: Solid
Participant2="solid-calculix"
cd ${Participant2}
    # Clean the case
    rm -fv *.log
    rm -fv tube.cvg
    rm -fv tube.dat
    rm -fv tube.frd
    rm -fv tube.sta
    rm -fv tube.12d
    rm -fv spooles.out
    rm -fv precice-*.log
    rm -fv precice-*-events.json
cd ..

# Remove the preCICE address files
rm -rfv precice-run

echo "Cleaning complete!"
#------------------------------------------------------------------------------
