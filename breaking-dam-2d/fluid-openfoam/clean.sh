#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_openfoam .
rm -fv 0/alpha.water
