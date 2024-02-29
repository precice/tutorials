#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh
rm -r preCICE-output

clean_calculix .
