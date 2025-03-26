#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
  echo "No input argument provided. Macro solver is launched in serial"
  ./macro_dumux params.input
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    ./macro_dumux params.input
    ;;
  p)
    mpiexec -n "$2" macro_dumux params.input
    ;;
  *)
    usage
    ;;
  esac
done

close_log
