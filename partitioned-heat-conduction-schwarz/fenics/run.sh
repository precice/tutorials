#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-d] [-n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	usage
fi

# Select appropriate case
while getopts ":dn" opt; do
  case ${opt} in
  d)
    python3 heat.py -d --error-tol 10e-3
    ;;
  n)
    python3 heat.py -n --error-tol 10e-3
    ;;
  *)
    usage
    ;;
  esac
done
