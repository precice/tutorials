#!/bin/sh
set -e -u

. ../../tools/log.sh

usage() { echo "Usage: cmd [-d] [-n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
        log usage
fi

log python3 -m venv .venv
log . .venv/bin/activate
log pip install -r requirements.txt

while getopts ":dn" opt; do
  case ${opt} in
  d)
    log rm -rf Dirichlet-*.vtk
    log NUTILS_RICHOUTPUT=no python3 heat.py --side=Dirichlet
    ;;
  n)
    log rm -rf Neumann-*.vtk
    log NUTILS_RICHOUTPUT=no python3 heat.py --side=Neumann
    ;;
  *)
    log usage
    ;;
  esac
done

close_log
