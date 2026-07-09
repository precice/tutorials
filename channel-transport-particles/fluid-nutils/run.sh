#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]; then
    if [ ! -d ".venv" ]; then
        python3 -m venv .venv
        source .venv/bin/activate
        pip install -r requirements.txt && pip freeze > pip-installed-packages.log
    else
        source .venv/bin/activate
    fi
fi

python3 fluid.py

# System tests: Keep the container and the respective network alive till the end.
if [[ -v PRECICE_TUTORIALS_TESTING ]]; then
    echo "Waiting for the Particles participant to finish..."
    if [ ! -f "../particles-participant-finished.log" ]; then
        inotifywait -e create,modify,attrib --include '/particles-participant-finished\.log$' -qq ..
    fi
fi

close_log
