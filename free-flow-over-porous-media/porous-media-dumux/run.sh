#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
  echo "No input argument provided. Porous medium solver is launched in serial"
  ./porous_media_dumux params.input
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    ./porous_media_dumux params.input
    ;;
  p)
    mpiexec -n "$2" porous_media_dumux params.input
    ;;
  *)
    usage
    ;;
  esac
done

close_log
