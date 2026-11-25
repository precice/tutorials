#!/usr/bin/env bash
set -e -u

python3 -m venv .venv
. .venv/bin/activate

pip install -r requirements.txt

python surrogate_workflow.py
