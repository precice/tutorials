#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-l] [-r]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	usage
fi

# Select appropriate case
while getopts ":lr" opt; do
  case ${opt} in
  l)
    python3 oscillator.py Mass-Left
    ;;
  r)
    python3 oscillator.py Mass-Right
    ;;
  *)
    usage
    ;;
  esac
done