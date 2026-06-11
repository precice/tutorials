#!/usr/bin/env bash
set -e -u

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    if [ ! -d .venv ]; then
        python3 -m venv --system-site-packages .venv
    fi
    . .venv/bin/activate
    pip install -r ../solver-fenicsx/requirements.txt
fi

python3 ../solver-fenicsx/heat.py Dirichlet --error-tol 10e-3
