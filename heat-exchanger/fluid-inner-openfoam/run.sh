#!/bin/sh
set -e -u

touch fluid-inner-openfoam.foam

../../tools/run-openfoam.sh "$@"
