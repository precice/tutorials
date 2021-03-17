#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

echo "Cleaning..."

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

# Participant 1: Fluid (OpenFOAM)
Participant1="fluid-openfoam"
cd ${Participant1}
    # Clean the case
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
    rm -fv \
        precice-*.log \
        precice-*-events.json
cd ..

# Participant 2: Solid1 (deal.II)
Participant2="solid-left-dealii"
cd  ./Solid1/dealii_output
    # Clean the case
    echo "Cleaning Solid1 case"
    rm -fv solution-*.vtk
cd ..
    rm -fv \
        precice-*.log \
        precice-*-events.json \
        solution-*.vtk
cd ..

# Participant 3: Solid2 (deal.II)
Participant3="solid-right-dealii"
cd  ./Solid2/dealii_output
    # Clean the case
    echo "Cleaning Solid2 case"
    rm -fv solution-*.vtk
cd ..
    rm -fv \
        precice-*.log \
        precice-*-events.json \
        solution-*.vtk
cd ..

# Remove the preCICE-related log files
echo "Deleting the preCICE log files..."
rm -fv \
    precice-*.log \
    precice-*-events.json

rm -rfv precice-run
rm -rfv precice-output


echo "Cleaning complete!"
#------------------------------------------------------------------------------
