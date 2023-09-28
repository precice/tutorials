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
    python3 heat.py -d
    ;;
  n)
    python3 heat.py -n
    ;;
  *)
    usage
    ;;
  esac
done
