#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ "${PRECICE_TUTORIALS_VENV:-true}" = true ]
then
    python3 -m venv .venv
    . .venv/bin/activate
    pip install -r requirements.txt && pip freeze > pip-installed-packages.log
fi

python3 transport.py "$@"

close_log
