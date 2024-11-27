#!/bin/sh
set -e -u

python3 -m venv .venv
. .venv/bin/activate
pip install ../solver-python

if [ $# -eq 0 ]; then
  growing A
else
  mpirun -n "$@" growing A
fi
