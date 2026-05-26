#!/usr/bin/env sh
set -e # Not setting -u as it gets triggered by the OpenFOAM RunFunctions

# Prepare an (intentionally empty) .foam file for the ParaView OpenFOAM reader
CASENAME="$(pwd | xargs basename)"
touch "$CASENAME.foam"

# Keep a backup of the files to modify
echo "backing up the original files (copies: 0/U.orig, system/controlDict.orig, constant/dynamicMeshDict.orig)"
cp 0/U 0/U.orig
cp system/controlDict system/controlDict.orig
cp constant/dynamicMeshDict constant/dynamicMeshDict.orig

# Modify code for foam-extend
echo "modifying everything now"
sed -i "s/noSlip;/noSlipWall;/g" 0/U
sed -i "s,application     pimpleFoam;,//application     pimpleFoam;,g" system/controlDict
sed -i "s,// application     pimpleDyMFoam;,application     pimpleDyMFoam;,g" system/controlDict
sed -i '41i\ \ \ \ "liblduSolvers.so"' system/controlDict
sed -i '41i\ \ \ \ "libforces.so"' system/controlDict
sed -i "s,writeCompression    off,writeCompression    uncompressed,g" system/controlDict
sed -i "s,\/\* uncomment,// FOAMEXTENDBEGIN,g" system/controlDict
sed -i "s,\*\/ // foam-extend,// FOAMEXTENDEND,g" system/controlDict

sed -i "s,\/\* uncomment,// FOAMEXTENDBEGIN,g" system/fvSchemes
sed -i "s,\*\/ // foam-extend,// FOAMEXTENDEND,g" system/fvSchemes

sed -i "s,\/\* uncomment,// FOAMEXTENDBEGIN,g" system/fvSolution
sed -i "s,\*\/ // foam-extend,// FOAMEXTENDEND,g" system/fvSolution

sed -i "s/libfvMotionSolvers\./libfvMotionSolver\./g" constant/dynamicMeshDict
sed -i "s,\/\* uncomment,// FOAMEXTENDBEGIN,g" constant/dynamicMeshDict
sed -i "s,\*\/ // foam-extend,// FOAMEXTENDEND,g" constant/dynamicMeshDict

# OpenFOAM run functions: getApplication, getNumberOfProcessors
# shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
. "${WM_PROJECT_DIR}/bin/tools/RunFunctions"
solver=$(getApplication | cut -f 1 -d " " | sed '\~//~d')
if [ "${1:-}" = "-parallel" ]; then
    procs=$(getNumberOfProcessors)
    decomposePar -force
    mpirun -np "${procs}" "${solver}" -parallel
    reconstructPar
else
    ${solver}
fi

# Reverse code for OpenFOAM
#rm -rf constant/polyMesh
#sed -i "s/noSlipWall;/noSlip;/g" 0/U
#sed -i "s,application     pimpleDyMFoam;,// application     pimpleDyMFoam;,g" system/controlDict
#sed -i "s,//application     pimpleFoam;,application     pimpleFoam;,g" system/controlDict
#sed -i '/   "liblduSolvers.so"/d' system/controlDict
#sed -i '/   "libforces.so/d' system/controlDict
#sed -i "s,writeCompression    uncompressed,writeCompression    off,g" system/controlDict
#sed -i "s,// FOAMEXTENDBEGIN,\/\* uncomment,g" system/controlDict
#sed -i "s,// FOAMEXTENDEND,\*\/ // foam-extend,g" system/controlDict
#
#sed -i "s,// FOAMEXTENDBEGIN,\/\* uncomment,g" system/fvSchemes
#sed -i "s,// FOAMEXTENDEND,\*\/ // foam-extend,g" system/fvSchemes
#
#sed -i "s,// FOAMEXTENDBEGIN,\/\* uncomment,g" system/fvSolution
#sed -i "s,// FOAMEXTENDEND,\*\/ // foam-extend,g" system/fvSolution
#
#sed -i "s/libfvMotionSolver\./libfvMotionSolvers\./g" constant/dynamicMeshDict
#sed -i "s,// FOAMEXTENDBEGIN,\/\* uncomment,g" constant/dynamicMeshDict
#sed -i "s,// FOAMEXTENDEND,\*\/ // foam-extend,g" constant/dynamicMeshDict
