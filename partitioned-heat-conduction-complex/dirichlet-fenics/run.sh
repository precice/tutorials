#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 ../solver-fenics/heat.py -d -i complex

close_log
