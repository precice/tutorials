#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-d] [-n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
        usage
fi

while getopts ":dn" opt; do
  case ${opt} in
  d)
    rm -rf Dirichlet-*.vtk
    ./heat D 10e-6
    ;;
  n)
    rm -rf Neumann-*.vtk
    ./heat N 10e-6
    ;;
  *)
    usage
    ;;
  esac
done
