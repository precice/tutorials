#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

usage() { echo "Usage: cmd [-s] [-d]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	usage
fi

# Select appropriate case
while getopts ":sd" opt; do
  case ${opt} in
  s)
    python3 volume-coupled-diffusion.py --source
    ;;
  d)
    python3 volume-coupled-diffusion.py --drain
    ;;
  *)
    usage
    ;;
  esac
done

close_log
