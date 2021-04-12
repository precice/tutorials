#!/bin/sh
set -e -u

if [ "${1:-}" = "-parallel" ]; then
    mpirun -n 2 SU2_CFD euler_config_coupled.cfg
else
    SU2_CFD euler_config_coupled.cfg
fi
