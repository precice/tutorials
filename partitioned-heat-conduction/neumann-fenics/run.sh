#!/usr/bin/env bash
set -e -u

python3 -m venv --system-site-packages ../solver-fenics/.venv
. ../solver-fenics/.venv/bin/activate
pip install -r ../solver-fenics/requirements.txt

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 ../solver-fenics/heat.py Neumann

close_log
