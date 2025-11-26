#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    python3 -m venv --system-site-packages .venv
    . .venv/bin/activate
    pip install -r ../solver-fenics/requirements.txt
fi

python3 ../solver-fenics/volume-coupled-diffusion.py --drain

close_log
