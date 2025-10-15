#!/usr/bin/env bash
set -e -u

if [ -d ".venv" ]; then
	echo "Using existing virtual environment"
	. .venv/bin/activate
else
	echo "Creating new virtual environment"
	python3 -m venv .venv
	. .venv/bin/activate
	pip install -r requirements.txt && pip freeze > pip-installed-packages.log
fi

if [ -f "../initial_condition.npz" ]; then
	:
else
	echo "Error, missing initial condition file."
	echo "Run 'python3 ../generate_ic.py --epoch <seed>' to create one."
	exit 1
fi

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

echo "[preCICE] Waiting for Dirichlet participant..."
python3 solver.py

close_log