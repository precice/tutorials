#!/bin/sh
set -e -u

. ../../tools/log.sh

while getopts ":dn" opt; do
  case ${opt} in
  d)
    log python3 heat.py -d --error-tol 10e-3
    ;;
  n)
    log python3 heat.py -n --error-tol 10e-3
    ;;
  \?)
    log echo "Usage: cmd [-d] [-n]"
    ;;
  esac
done

close_log
