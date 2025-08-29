#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 -m venv .venv
. .venv/bin/activate
pip install ../solver-python

if [ $# -eq 0 ]; then
  growing B
else
  mpirun -n "$@" growing B
fi

close_log