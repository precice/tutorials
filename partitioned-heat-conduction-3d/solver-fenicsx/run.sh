#!/usr/bin/env bash
set -e -u

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    if [ ! -d ".venv" ]; then
        python3 -m venv --system-site-packages .venv
        . .venv/bin/activate
        pip install ../../..
        pip install -r .requirements.txt
    else
        . .venv/bin/activate
    fi
fi

while getopts ":dn" opt; do
  case ${opt} in
  d)
    python3 heat.py Dirichlet --error-tol 10e-3
    ;;
  n)
    python3 heat.py Neumann --error-tol 10e-3
    ;;
  \?)
    echo "Usage: cmd [-d] [-n]"
    ;;
  esac
done
