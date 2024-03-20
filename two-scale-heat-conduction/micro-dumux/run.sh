#!/bin/sh
set -e -u

. ../../tools/log.sh

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
  log echo "No input argument provided. Micro Manager is launched in serial"
  log python3 run_micro_manager.py params.input
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    log python3 run_micro_manager.py params.input
    ;;
  p)
    log mpiexec -n "$2" python3 run_micro_manager.py params.input
    ;;
  *)
    log usage
    ;;
  esac
done

close_log
