#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh
rm -r preCICE-CCX-output

clean_calculix .
