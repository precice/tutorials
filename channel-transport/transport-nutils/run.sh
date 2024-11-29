#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python3 transport.py "$@"

close_log
