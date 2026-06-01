#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    if [ ! -d .venv ]; then
        python3 -m venv --system-site-packages .venv
    fi
    . .venv/bin/activate
    pip install -r requirements.txt && pip freeze > pip-installed-packages.log
fi

SU2_preCICE_CHT.py -f laminar_config_unsteady.cfg -r --parallel

close_log
