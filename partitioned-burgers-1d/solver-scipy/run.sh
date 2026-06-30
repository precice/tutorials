#!/usr/bin/env bash
set -e -u

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    if [ ! -d ".venv" ]; then
        python3 -m venv .venv
        source .venv/bin/activate
        pip install -r requirements.txt && pip freeze > pip-installed-packages.log
    else
        source .venv/bin/activate
    fi
fi

if [ ! -f "../initial_condition.npz" ]; then
	echo "Generating initial condition..."
	python3 ../utils/generate_ic.py
fi

# Run the monolithic solver
# The 'None' argument tells the solver to run monolithic (preCICE participant name is none, run without preCICE)
# Append any additional arguments that this script has been called with.
python3 solver.py None "$@"
