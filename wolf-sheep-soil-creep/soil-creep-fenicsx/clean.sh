#!/usr/bin/env sh
set -e -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

rm -rfv ./output/

. ../../tools/cleaning-tools.sh

clean_fenicsx .
