#!/bin/sh
set -e -u

usage() { echo "Usage: cmd [-l] [-r]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
	usage
fi

# Select appropriate case
while getopts ":lr" opt; do
  case ${opt} in
  l)
    python3 runner.py ./MassLeft/fmi-settings.json MassLeft/precice-settings.json
    python3 calculate-error.py MassLeft/fmi-settings.json MassLeft/precice-settings.json MassRight/fmi-settings.json MassRight/precice-settings.json Mass-Left

    ;;
  r)
    python3 runner.py MassRight/fmi-settings.json MassRight/precice-settings.json
    python3 calculate-error.py MassLeft/fmi-settings.json MassLeft/precice-settings.json MassRight/fmi-settings.json MassRight/precice-settings.json Mass-Right

    ;;
  *)
    usage
    ;;
  esac
done
