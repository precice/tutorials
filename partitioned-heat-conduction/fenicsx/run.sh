#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

while getopts ":dn" opt; do
  case ${opt} in
  d)
    python3 heat.py -d --error-tol 10e-3
    ;;
  n)
    python3 heat.py -n --error-tol 10e-3
    ;;
  \?)
    echo "Usage: cmd [-d] [-n]"
    ;;
  esac
done

close_log
