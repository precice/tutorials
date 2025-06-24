#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r ../solver-fenics/requirements.txt

python3 ../solver-fenics/heat.py Right

close_log
