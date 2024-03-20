#!/bin/sh
set -e -u

. ../../tools/log.sh

usage() { echo "Usage: run.sh [-modal]" 1>&2; exit 1; }

# There must be either 0 arguments or 1, which is modal.
# Send an error otherwise
if [ $# -ge 2 ] || { [ $# -eq 1 ] && [ "$1" != "-modal" ]; }; then
    log usage
fi

# No arg => regular simulation. Otherwise, it's modal
if [ $# -eq 0 ]; then
    log ccx_preCICE -i flap -precice-participant Solid
else
    log ccx_preCICE -i frequency
    mv frequency.eig flap_modal.eig
    log ccx_preCICE -i flap_modal -precice-participant Solid
fi

close_log
