#!/bin/sh
set -e -u

. ../../tools/log.sh

if [ "${1:-}" = "-parallel" ]; then
    log mpirun -n 2 SU2_CFD euler_config_coupled.cfg
else
    log SU2_CFD euler_config_coupled.cfg
fi

close_log
