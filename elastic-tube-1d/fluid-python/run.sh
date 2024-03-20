#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 ./FluidSolver.py ../precice-config.xml

close_log
