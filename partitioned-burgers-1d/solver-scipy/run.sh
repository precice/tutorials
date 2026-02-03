#!/usr/bin/env bash
set -e -u

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    python3 -m venv .venv
    . .venv/bin/activate
    pip install -r requirements.txt && pip freeze > pip-installed-packages.log
fi

if [ ! -f "../initial_condition.npz" ]; then
	echo "Generating initial condition..."
	python3 ../utils/generate_ic.py
fi

python3 solver.py None # run monolithic reference solution
