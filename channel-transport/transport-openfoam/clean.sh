#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_openfoam .
rm -rf ./0/
rm -rf ./Make/linux*
