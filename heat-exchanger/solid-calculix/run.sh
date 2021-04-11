#!/bin/sh
set -e -u

if [ ! -f all.msh ]; then
    echo "Mesh files not found. Use the Download_meshes script to download them."
    exit
fi

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
ccx_preCICE -i solid -precice-participant Solid
