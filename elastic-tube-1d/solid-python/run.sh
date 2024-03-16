#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 ./SolidSolver.py ../precice-config.xml

close_log
