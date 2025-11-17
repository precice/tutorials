#!/usr/bin/env sh
set -e -u

if [ "${PRECICE_TUTORIALS_VENV:-true}" = true ]
then
    python3 -m venv .venv
    . .venv/bin/activate
    pip install -r requirements.txt && pip freeze > pip-installed-packages.log
fi
python3 fake.py
