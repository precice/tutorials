#!/bin/sh
set -e -u

. ../../tools/log.sh

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
log ccx_preCICE -i tube -precice-participant Solid

close_log
