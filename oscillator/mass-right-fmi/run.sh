#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -f Oscillator.fmu ]; then
  echo "Building a copy of ../solver-fmi/fmu..."
  (
    cp -r ../solver-fmi _solver-fmi-copy
    cd _solver-fmi-copy/fmu
    rm -rf build && mkdir build && cd build
    # Both FMI_VERSION=3 and FMI_VERSION=2 are supported
    cmake -DFMI_TYPE=CS -DFMI_VERSION=3 ..
    make
    cp ./Oscillator.fmu ../../..
  )
  rm -r _solver-fmi-copy
fi

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]; then
    if [ ! -d ".venv" ]; then
        python3 -m venv .venv
        source .venv/bin/activate
        pip install -r requirements.txt && pip freeze > pip-installed-packages.log
    else
        source .venv/bin/activate
    fi
fi

fmiprecice fmi-settings.json precice-settings.json
python3 ../solver-fmi/calculate-error.py ../mass-left-fmi/fmi-settings.json ../mass-left-fmi/precice-settings.json ../mass-right-fmi/fmi-settings.json ../mass-right-fmi/precice-settings.json Mass-Right

close_log
