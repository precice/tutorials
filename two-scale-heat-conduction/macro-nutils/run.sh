#!/bin/sh
set -e -u

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
NUTILS_RICHOUTPUT=no python3 macro.py
