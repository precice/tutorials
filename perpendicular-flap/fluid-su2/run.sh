#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ $# -eq 0 ]
then
    echo "Installing dependencies in a Python virtual environment"
    python3 -m venv --system-site-packages .venv
    . .venv/bin/activate
    pip install -r requirements.txt && pip freeze > pip-installed-packages.log
fi

SU2_preCICE_FSI.py -f euler_config_unsteady.cfg --parallel

close_log
