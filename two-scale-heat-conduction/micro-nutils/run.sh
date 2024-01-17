#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt

# Check if no input argument was provided
if [ -z "$*" ] ; then
        echo "No input argument provided. Micro Manager is launched in serial"
	python3 run-micro-problems.py
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    python3 run-micro-problems.py
    ;;
  p)
    mpiexec -n "$2" python3 run-micro-problems.py
    ;;
  *)
    usage
    ;;
  esac
done
