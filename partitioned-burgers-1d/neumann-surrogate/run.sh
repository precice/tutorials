#!/usr/bin/env bash
set -e -u

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt

if [ ! -f "../initial_condition.npz" ]; then
	echo "Generating initial condition..."
	python3 ../utils/generate_ic.py
fi

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 solver.py

close_log