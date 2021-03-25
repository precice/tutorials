#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_su2 .
rm -fv fluid-su2.log
