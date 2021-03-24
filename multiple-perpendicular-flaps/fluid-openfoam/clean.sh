#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_openfoam .
rm -rfv ./0/ # here, 0.orig/ is copied to 0/