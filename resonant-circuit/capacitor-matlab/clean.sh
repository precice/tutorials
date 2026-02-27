#!/bin/sh
set -e -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

rm ./*.png ./*.mat

. ../../tools/cleaning-tools.sh

clean_matlab .
