#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -rfv ./output/

clean_precice_logs .
clean_case_logs .
