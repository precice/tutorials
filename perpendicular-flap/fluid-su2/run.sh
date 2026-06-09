#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]; then
    if [ ! -d ".venv" ]; then
        python3 -m venv --system-site-packages .venv
        source .venv/bin/activate
        pip install -r requirements.txt && pip freeze > pip-installed-packages.log
    else
        source .venv/bin/activate
    fi
fi

SU2_preCICE_FSI.py -f euler_config_unsteady.cfg --parallel

close_log
