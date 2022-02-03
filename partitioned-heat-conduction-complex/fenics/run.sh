#!/bin/sh
set -e -u

while getopts ":dn" opt; do
  case ${opt} in
  d)
    python3 heat.py -d -a
    ;;
  n)
    python3 heat.py -n -a
    ;;
  \?)
    echo "Usage: cmd [-d] [-n]"
    ;;
  esac
done
