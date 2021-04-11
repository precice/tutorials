#!/bin/sh
set -e -u

touch fluid-outer-openfoam.foam

../../tools/run-openfoam.sh "$@"
