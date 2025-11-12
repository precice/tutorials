#!/usr/bin/env sh
set -e -u

if [ $# -eq 0 ] 
then
    echo "Installing dependencies in a Python virtual environment"
    python3 -m venv .venv
    . .venv/bin/activate
    pip install -r requirements.txt && pip freeze > pip-installed-packages.log
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
python3 fake.py
