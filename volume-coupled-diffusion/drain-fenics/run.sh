#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ $# -eq 0 ]
then
    echo "Installing dependencies in a Python virtual environment"
    python3 -m venv --system-site-packages .venv
    . .venv/bin/activate
    pip install -r ../solver-fenics/requirements.txt
fi

python3 ../solver-fenics/volume-coupled-diffusion.py --drain

close_log
