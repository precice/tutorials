#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    if [ ! -d .venv ]; then
        python3 -m venv .venv
    fi
    . .venv/bin/activate
    pip install -r requirements.txt && pip freeze > pip-installed-packages.log
fi

python3 source.py

# System tests: Keep the container and the respective network alive till the end.
echo "Waiting for the Fluid participant to finish..."
inotifywait -e create,attrib -qq ../fluid-participant-finished.log

close_log
