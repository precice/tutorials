#!/usr/bin/env bash
set -e -u

cd "$(cd "$(dirname "$0")" && pwd)"
. ../../tools/cleaning-tools.sh

clean_precice_logs .
clean_case_logs .
