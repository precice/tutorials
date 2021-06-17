#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-a] [-i]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	usage
fi

# Select appropriate case
while getopts ":ai" opt; do
  case ${opt} in
  a)
    python3 gierer-meinhardt.py --activator
    ;;
  i)
    python3 gierer-meinhardt.py --inhibitor
    ;;
  *)
    usage
    ;;
  esac
done
