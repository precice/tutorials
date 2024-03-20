#!/bin/bash
set -e -u

. ../../tools/log.sh

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
ccx_preCICE -i beam2 -precice-participant Calculix2

close_log
