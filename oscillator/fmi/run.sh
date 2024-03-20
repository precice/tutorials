#!/bin/sh
set -e -u

. ../../tools/log.sh

usage() { echo "Usage: cmd [-l] [-r]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	log usage
fi

if [ ! -f Oscillator.fmu ]; then
  log cd fmu
  log rm -rf build
  log mkdir build
  log cd build
  log # Both FMI_VERSION=3 and FMI_VERSION=2 are supported
  log cmake -DFMI_TYPE=CS -DFMI_VERSION=3 ..
  log make
  log cp ./Oscillator.fmu ../..
  log cd ../../
fi

# Select appropriate case
while getopts ":lr" opt; do
  case ${opt} in
  l)
    log fmiprecice ./MassLeft/fmi-settings.json MassLeft/precice-settings.json
    log python3 calculate-error.py MassLeft/fmi-settings.json MassLeft/precice-settings.json MassRight/fmi-settings.json MassRight/precice-settings.json Mass-Left

    ;;
  r)
    log fmiprecice MassRight/fmi-settings.json MassRight/precice-settings.json
    log python3 calculate-error.py MassLeft/fmi-settings.json MassLeft/precice-settings.json MassRight/fmi-settings.json MassRight/precice-settings.json Mass-Right

    ;;
  *)
    log usage
    ;;
  esac
done

close_log
