#!/bin/sh
set -e -u

. ../../tools/log.sh

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

log python3 -m venv .venv
log . .venv/bin/activate
log pip install -r requirements.txt

# Check if no input argument was provided
if [ -z "$*" ] ; then
  log echo "No input argument provided. Micro Manager is launched in serial"
  log python3 run-micro-problems.py
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    log python3 run-micro-problems.py
    ;;
  p)
    log mpiexec -n "$2" python3 run-micro-problems.py
    ;;
  *)
    log usage
    ;;
  esac
done

close_log
