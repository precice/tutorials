#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
  echo "No input argument provided. Free flow solver is launched in serial"
  ./free_flow_dumux params.input
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    ./free_flow_dumux params.input
    close_log
    ;;
  p)
    echo "Participant is not started, only serial execution is available."
    ;;
  *)
    usage
    ;;
  esac
done