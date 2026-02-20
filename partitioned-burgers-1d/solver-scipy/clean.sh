#!/usr/bin/env sh
set -e -u

# shellcheck disable=SC1091
. ../../tools/cleaning-tools.sh

clean_case_logs .
rm -f full_domain.npz
