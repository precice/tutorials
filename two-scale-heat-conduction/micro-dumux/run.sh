#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
        echo "No input argument provided. Micro Manager is launched in serial"
        cd build-cmake/appl/
	python3 run_micro_manager.py params.input
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    cd build-cmake/appl/
    python3 run_micro_manager.py params.input
    ;;
  p)
    cd build-cmake/appl/
    mpiexec -n "$2" python3 run_micro_manager.py params.input
    ;;
  *)
    usage
    ;;
  esac
done
