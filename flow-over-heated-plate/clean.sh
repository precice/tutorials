#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

echo "Cleaning..."

clean_case(){
   cd ${1}
    # Clean the result and auxiliary files
    rm -fv *.vtk
    rm -fv ${1}.log
    rm -fv precice-*.log \
    rm -fv precice-postProcessingInfo.log \
    rm -fv precice-*-events.json
}

finished_clean(){
cd ..
}

# Nutils
Name="solid-nutils"
clean_case ${Name}
finished_clean

# Openfoam
Name="fluid-openfoam"
clean_case ${Name}
    # Clean specialized files
    if [ -n "${WM_PROJECT}" ]; 
    then
        . $WM_PROJECT_DIR/bin/tools/CleanFunctions
        cleanCase
        rm -rfv ./0
        touch ${Name}.foam
    fi
finished_clean

# Openfoam
Name="solid-openfoam"
clean_case ${Name}
    # Clean specialized files
    if [ -n "${WM_PROJECT}" ];
    then
        . $WM_PROJECT_DIR/bin/tools/CleanFunctions
        cleanCase
        rm -rfv ./0
        touch ${Name}.foam
    fi
finished_clean

# Fenics
Name="solid-fenics"
clean_case ${Name}
    # Clean specialized files
    rm -fv *.pvd
    rm -fv FSI-S/*
    rm -fv spooles.out
finished_clean


# Remove the preCICE address files
rm -rfv precice-run
rm -fv .*.address

echo "Cleaning complete!"
#------------------------------------------------------------------------------
