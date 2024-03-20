#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 ../solver-fenics/heat.py Right

close_log
