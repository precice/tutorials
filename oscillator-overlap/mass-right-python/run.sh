#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 ../solver-python/oscillator.py Mass-Right

close_log
