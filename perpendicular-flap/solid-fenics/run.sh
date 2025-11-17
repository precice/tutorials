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

python3 solid.py

close_log
