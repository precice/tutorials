#!/bin/sh
set -e -u

usage() { echo "Usage: run.sh [-modal]" 1>&2; exit 1; }

if [ $# -ge 2 ] || [ $# -eq 1 -a $1 != "-modal" ]; then
    usage
fi
#ccx_preCICE -i flap -precice-participant Solid


if [ ${1-nonmodal} = "-modal" ]; then
    echo "Modal"
    echo ${1-nonmodal}
else
    echo "NonModal"
    echo ${1-nonmodal}
fi