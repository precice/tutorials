#!/bin/sh
set -e -u

# Set the working directory to be the one where this script is located
cd "$(dirname "$0")"

python3 launch_unsteady_CHT_FlatPlate.py -f laminar_config_unsteady.cfg --parallel