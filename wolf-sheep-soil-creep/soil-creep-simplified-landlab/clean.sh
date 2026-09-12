#!/usr/bin/env sh
set -e -u

rm -rfv ./output/

. ../../tools/cleaning-tools.sh

clean_precice_logs .
clean_case_logs .
