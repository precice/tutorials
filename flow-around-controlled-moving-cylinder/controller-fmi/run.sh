#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -f PIDcontroller.fmu ]; then
  cd fmu
  rm -rf build
  mkdir build
  cd build
  cmake -DFMI_TYPE=CS -DFMI_VERSION=3 ..
  make
  cp ./PIDcontroller.fmu ../..
  cd ../../
fi

python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt && pip freeze > pip-installed-packages.log

fmiprecice ./fmi-settings.json ./precice-settings.json

close_log
