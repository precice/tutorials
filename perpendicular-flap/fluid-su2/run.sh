#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ "${1:-}" = "-parallel" ]; then
    mpirun -n 2 SU2_CFD euler_config_coupled.cfg
else
    SU2_CFD euler_config_coupled.cfg
fi

close_log
