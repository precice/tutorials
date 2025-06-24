#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt && pip freeze > pip-installed-packages.log

rm -rf Dirichlet-*.vtk
NUTILS_RICHOUTPUT=no python3 heat.py

close_log
