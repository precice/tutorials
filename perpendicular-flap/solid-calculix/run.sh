#!/bin/sh
set -e -u

usage() { echo "Usage: run.sh [-modal]" 1>&2; exit 1; }

if [ $# -ge 2 ] || [ $# -eq 1 -a $1 != "-modal" ]; then
    usage
fi

if [ ${1-nonmodal} = "-modal" ]; then
    ccx_preCICE -i frequency
    mv frequency.eig flap_modal.eig
    ccx_preCICE -i flap_modal -precice-participant Solid
else
    ccx_preCICE -i flap -precice-participant Solid
fi