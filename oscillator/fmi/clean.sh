#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -rfv ./output/
rm -rfv ../precice-run/

clean_precice_logs .
