#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
ccx_preCICE -i tube -precice-participant Solid

close_log
