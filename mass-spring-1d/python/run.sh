#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-m] [-n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	usage
fi

# Select appropriate case
while getopts ":mn" opt; do
  case ${opt} in
  m)
    python3 mass-spring.py MassOne
    ;;
  n)
    python3 mass-spring.py MassTwo
    ;;
  *)
    usage
    ;;
  esac
done