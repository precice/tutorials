#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 ../solver-python/oscillator.py Mass-Left

close_log
