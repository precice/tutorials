#!/usr/bin/env bash
set -e -u

solver_path="../solver-scipy"

if [ -d "$solver_path/.venv" ]; then
	echo "Using existing virtual environment"
	. "$solver_path/.venv/bin/activate"
else
	echo "Creating new virtual environment"
	python3 -m venv "$solver_path/.venv"
	. "$solver_path/.venv/bin/activate"
	pip install -r $solver_path/requirements.txt && pip freeze > $solver_path/pip-installed-packages.log
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
python3 "$solver_path/solver.py" Neumann

close_log
