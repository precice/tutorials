#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    if [ ! -d ".venv" ]; then
        python3 -m venv --system-site-packages .venv
        source .venv/bin/activate
        pip install -r ../solver-fenics/requirements.txt
    else
        source .venv/bin/activate
    fi
fi

python3 ../solver-fenics/volume-coupled-diffusion.py --drain

close_log
