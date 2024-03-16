#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 ./FluidSolver.py ../precice-config.xml

close_log
