#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
        usage
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    python3 run-micro-problems.py
    ;;
  p)
    mpiexec -n $2 python3 run-micro-problems.py
    ;;
  *)
    usage
    ;;
  esac
done
