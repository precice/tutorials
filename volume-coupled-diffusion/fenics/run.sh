#!/bin/sh
set -e -u

. ../../tools/log.sh

usage() { echo "Usage: cmd [-s] [-d]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	log usage
fi

# Select appropriate case
while getopts ":sd" opt; do
  case ${opt} in
  s)
    log python3 volume-coupled-diffusion.py --source
    ;;
  d)
    log python3 volume-coupled-diffusion.py --drain
    ;;
  *)
    log usage
    ;;
  esac
done

close_log
