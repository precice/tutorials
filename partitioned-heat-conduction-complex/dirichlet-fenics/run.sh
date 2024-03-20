#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 ../solver-fenics/heat.py -d -i complex

close_log
