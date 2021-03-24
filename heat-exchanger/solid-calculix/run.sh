#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

# This script prepares and runs all the participants in one terminal,
# forwarding the solvers' output to log files.
# You need to first download the mesh files separately using the Download_meshes script.
# The script "Allclean" cleans-up the result and log files.

if [ ! -f all.msh ]; then
    echo "Mesh files not found. Use the Download_meshes script to download them."
    exit
fi

Participant="Solid"
Solver="ccx_preCICE"

# Run and get the process id
echo "Starting the ${Participant} participant..."
export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
${Solver} -i solid -precice-participant ${Participant}
