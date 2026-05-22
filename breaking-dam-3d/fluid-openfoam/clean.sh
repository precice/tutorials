#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh
rm -r 0

clean_openfoam .
