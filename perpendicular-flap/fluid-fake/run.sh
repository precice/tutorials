#!/usr/bin/env bash
set -e -u

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    python3 -m venv .venv
    . .venv/bin/activate
    pip install -r requirements.txt && pip freeze > pip-installed-packages.log
fi
python3 fake.py
