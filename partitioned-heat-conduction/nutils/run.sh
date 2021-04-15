#!/bin/sh
set -e -u

while getopts ":dn" opt; do
  case ${opt} in
  d)
    rm -rf Dirichlet-*.vtk
    python3 heat.py --side=Dirichlet
    ;;
  n)
    rm -rf Neumann-*.vtk
    python3 heat.py --side=Neumann
    ;;
  \?)
    echo "Usage: cmd [-d] [-n]"
    ;;
  esac
done
