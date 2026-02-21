#!/usr/bin/env sh
set -e -u

cd "$(cd "$(dirname "$0")" && pwd)"
. ../../tools/cleaning-tools.sh

clean_openfoam .
