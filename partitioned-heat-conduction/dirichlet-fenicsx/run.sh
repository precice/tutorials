#!/bin/sh
set -e -u

if [ ! -d .venv ]; then
        python3 -m venv --system-site-packages .venv
    fi
. .venv/bin/activate
pip install -r ../solver-fenicsx/requirements.txt

python3 ../solver-fenicsx/heat.py Dirichlet
