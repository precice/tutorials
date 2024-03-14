#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-l] [-r]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	usage
fi

if [ ! -f Oscillator.fmu ]; then
  cd fmu
  rm -rf build
  mkdir build
  cd build
  # Both FMI_VERSION=3 and FMI_VERSION=2 are supported
  cmake -DFMI_TYPE=CS -DFMI_VERSION=2 ..
  make
  cp ./Oscillator.fmu ../..
  cd ../../
fi

# Select appropriate case
while getopts ":lr" opt; do
  case ${opt} in
  l)
    fmiprecice ./MassLeft/fmi-settings.json MassLeft/precice-settings.json
    python3 calculate-error.py MassLeft/fmi-settings.json MassLeft/precice-settings.json MassRight/fmi-settings.json MassRight/precice-settings.json Mass-Left

    ;;
  r)
    fmiprecice MassRight/fmi-settings.json MassRight/precice-settings.json
    python3 calculate-error.py MassLeft/fmi-settings.json MassLeft/precice-settings.json MassRight/fmi-settings.json MassRight/precice-settings.json Mass-Right

    ;;
  *)
    usage
    ;;
  esac
done
