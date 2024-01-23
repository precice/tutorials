#!/bin/sh
set -e -u

while getopts ":dn" opt; do
  case ${opt} in
  d)
    python3 heat.py -d -i complex
    ;;
  n)
    python3 heat.py -n -i complex
    ;;
  \?)
    echo "Usage: cmd [-d] [-n]"
    ;;
  esac
done
