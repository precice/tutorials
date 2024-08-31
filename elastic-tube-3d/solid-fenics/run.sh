#!/usr/bin/env bash
set -e -u

python3 -m venv --system-site-package .venv
. .venv/bin/activate
pip install -r requirements.txt

python3 solid.py
