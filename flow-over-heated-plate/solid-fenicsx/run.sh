#!/bin/sh
set -e -u

python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r requirements.txt

python3 solid.py
