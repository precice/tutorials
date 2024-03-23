#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

SU2_preCICE_FSI.py -f euler_config_unsteady.cfg --parallel

close_log
