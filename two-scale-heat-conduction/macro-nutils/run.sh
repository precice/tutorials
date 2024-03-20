#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 -m venv .venv
log . .venv/bin/activate
log pip install -r requirements.txt
log python3 macro.py richoutput=no

close_log
