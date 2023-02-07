#!/bin/sh
set -e -u

python3 launch_unsteady_FSI.py -f euler_config_unsteady.cfg --parallel

