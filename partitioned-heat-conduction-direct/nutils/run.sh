#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-d] [-n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
        usage
fi

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt

while getopts ":dn" opt; do
  case ${opt} in
  d)
    rm -rf Dirichlet-*.vtk
    NUTILS_RICHOUTPUT=no python3 heat.py --side=Dirichlet
    ;;
  n)
    rm -rf Neumann-*.vtk
    NUTILS_RICHOUTPUT=no python3 heat.py --side=Neumann
    ;;
  *)
    usage
    ;;
  esac
done
