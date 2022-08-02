#!/bin/bash

# Currently, the case has only been checked with OpenFOAM-v2012
if [[ "${WM_PROJECT}" != "OpenFOAM" || "${WM_PROJECT_VERSION}" != "v2012" ]]
then
    echo; echo "This case currently only runs in OpenFOAM-v2012"; echo
    exit 0
fi

blockMesh
touch solid-openfoam.foam

# Compile boundary condition
(cd solidDisplacementFoamForce && wmake libso)

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
