#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_calculix .
rm -fv Solid.log
