#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 -m venv .venv
. .venv/bin/activate
pip install ../solver-python
pip freeze ../solver-python/ > ../solver-python/pip-installed-packages.log

if [ $# -eq 0 ]; then
  growing A
else
  mpirun -n "$@" growing A
fi

close_log