#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r requirements.txt

SU2_preCICE_FSI.py -f euler_config_unsteady.cfg --parallel

close_log
