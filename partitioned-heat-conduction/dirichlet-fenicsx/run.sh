#!/bin/sh
set -e -u

python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r ../solver-fenicsx/requirements.txt

python3 ../solver-fenicsx/heat.py Dirichlet
