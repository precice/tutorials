#!/usr/bin/env bash
set -e -u

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 ../solver-python/oscillator.py Mass-Left

close_log
