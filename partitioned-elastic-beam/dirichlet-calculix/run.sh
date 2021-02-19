#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

# This script prepares and runs the CalculiX solver

# =============== Participant: Solid ===========================
Solver="ccx_preCICE"

    # Run and get the process id
    echo "Starting the CalculiX..."
    export OMP_NUM_THREADS=1
    export CCX_NPROC_EQUATION_SOLVER=1
    ${Solver} -i beam1 -precice-participant Calculix1 
