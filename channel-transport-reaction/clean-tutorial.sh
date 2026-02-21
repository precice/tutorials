#!/usr/bin/env sh
set -e -u

cd "$(cd "$(dirname "$0")" && pwd)"

# shellcheck disable=SC1091
. ../tools/cleaning-tools.sh

clean_tutorial .
clean_precice_logs .
rm -fv ./*.log
rm -fv ./*.vtu

