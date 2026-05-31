#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_openfoam .

# Necessary signal file for the system tests,
# to keep the Source container running till the end
rm -f ../fluid-participant-finished.log