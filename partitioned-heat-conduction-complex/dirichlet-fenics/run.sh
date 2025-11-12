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
else
    case "$1" in
        -[s]|--skip-setup)
            echo "Skipping setup: Assuming an already prepared Python environment."
            ;;
        *)
            echo "Usage: $0 [-s|--skip-setup]"
            ;;
    esac
fi

python3 ../solver-fenics/heat.py -d -i complex

close_log
