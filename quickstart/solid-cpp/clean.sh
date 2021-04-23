#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -rfv coupling-meshes
clean_precice_logs .
