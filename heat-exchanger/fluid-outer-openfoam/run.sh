#!/bin/sh
set -e -u

rm -rfv 0/ && cp -r 0.orig/ 0/
touch fluid-outer-openfoam.foam

../../tools/run-openfoam.sh "$@"