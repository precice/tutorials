#!/bin/bash

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Currently only checked with OpenFOAM-v2012
if [[ "${WM_PROJECT}" != "OpenFOAM" || "${WM_PROJECT_VERSION}" != "v2012" ]]
then
    echo; echo "This case currently only runs in OpenFOAM-v2012"; echo
    exit 0
fi
# Source solids4Foam scripts
source solids4FoamScripts.sh

# Create mesh
runApplication blockMesh

# Run solver
if [ "${1:-}" = "-parallel" ]; then
    procs=$(getNumberOfProcessors)
    runApplication decomposePar -force
    mpirun -np "${procs}" solids4Foam -parallel
    runApplication reconstructPar
else
    runApplication solids4Foam
fi

# Remove empty time directories created when using preCICE
solids4Foam::removeEmptyDirs
