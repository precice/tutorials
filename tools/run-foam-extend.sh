#!/usr/bin/env sh
set -e # Not setting -u as it gets triggered by the OpenFOAM RunFunctions

# Prepare an (intentionally empty) .foam file for the ParaView OpenFOAM reader
CASENAME="$(pwd | xargs basename)"
touch "$CASENAME.foam"

# Modify code for foam-extend
mkdir constant/polyMesh
cp system/blockMeshDict constant/polyMesh 
sed -i "s/noSlip;/noSlipWall;/g" 0/U
sed -i "s,application     pimpleFoam;,//application     pimpleFoam;,g" system/controlDict
sed -i "s,// application     pimpleDyMFoam;,application     pimpleDyMFoam;,g" system/controlDict
sed -i '41i\ \ \ \ "liblduSolvers.so"' system/controlDict
sed -i '41i\ \ \ \ "libforces.so"' system/controlDict
sed -i "s,writeCompression    off,writeCompression    uncompressed,g" system/controlDict

sed -i "s/libfvMotionSolvers\./libfvMotionSolver\./g" constant/dynamicMeshDict

# OpenFOAM run functions: getApplication, getNumberOfProcessors
# shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
. "${WM_PROJECT_DIR}/bin/tools/RunFunctions"
solver=$(getApplication)
if [ "${1:-}" = "-parallel" ]; then
    procs=$(getNumberOfProcessors)
    decomposePar -force
    mpirun -np "${procs}" "${solver}" -parallel
    reconstructPar
else
    ${solver}
fi

# Reverse code for OpenFOAM
rm -rf constant/polyMesh
sed -i "s/noSlipWall;/noSlip;/g" 0/U
sed -i "s,application     pimpleDyMFoam;,// application     pimpleDyMFoam;,g" system/controlDict
sed -i "s,//application     pimpleFoam;,application     pimpleFoam;,g" system/controlDict
sed -i '/   "liblduSolvers.so"/d' system/controlDict
sed -i '/   "libforces.so/d' system/controlDict
sed -i "s,writeCompression    uncompressed,writeCompression    off,g" system/controlDict

sed -i "s/libfvMotionSolver\./libfvMotionSolvers\./g" constant/dynamicMeshDict
