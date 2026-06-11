#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    if [ ! -d ".venv" ]; then
        python3 -m venv .venv
        source .venv/bin/activate
        pip install -r ../solver-scipy/requirements.txt && pip freeze > pip-installed-packages.log
    else
        source .venv/bin/activate
    fi
fi

if [ ! -f "../initial_condition.npz" ]; then
	echo "Generating initial condition..."
	python3 ../utils/generate_ic.py
fi

python3 ../solver-scipy/solver.py Neumann

close_log
