#!/usr/bin/env sh
set -e -u

cd "$(cd "$(dirname "$0")" && pwd)"
. ../../tools/cleaning-tools.sh

# since we work with a 0_orig folder here
rm -rf 0

clean_openfoam .
