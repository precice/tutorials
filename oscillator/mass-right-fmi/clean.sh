#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -rfv ./output/

clean_precice_logs .
