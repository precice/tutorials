#!/bin/sh
set -e -u

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 [fans|nasmat]"
    exit 1
fi

if [ "$1" = "fans" ]; then
    echo "Using the CalculiX adapter configuration file for the coupling with the micro-mechanics solver FANS"
    cp config_for_fans.yml config.yml
elif [ "$1" = "nasmat" ]; then
    echo "Using the CalculiX adapter configuration file for the coupling with the micro-mechanics solver NASMAT"
    cp config_for_nasmat.yml config.yml
else
    echo "Invalid option: $1. Use 'fans' or 'nasmat'."
    exit 1
fi

$HOME/calculix-adapter/bin/ccx_preCICE -i notch -precice-participant Meso-structure

rm config.yml
