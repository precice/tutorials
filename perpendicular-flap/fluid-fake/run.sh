#!/usr/bin/env sh
set -e -u

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python3 fake.py
