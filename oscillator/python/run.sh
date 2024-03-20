#!/bin/sh
set -e -u

. ../../tools/log.sh

usage() { echo "Usage: cmd [-l] [-r]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	log usage
fi

# Select appropriate case
while getopts ":lr" opt; do
  case ${opt} in
  l)
    log python3 oscillator.py Mass-Left
    ;;
  r)
    log python3 oscillator.py Mass-Right
    ;;
  *)
    log usage
    ;;
  esac
done

close_log
