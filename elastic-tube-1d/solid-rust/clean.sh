#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_precice_logs .
clean_case_logs .
