#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -f ../solver-fmi/Oscillator.fmu ]; then
  cd ../solver-fmi/fmu
  rm -rf build
  mkdir build
  cd build
  # Both FMI_VERSION=3 and FMI_VERSION=2 are supported
  cmake -DFMI_TYPE=CS -DFMI_VERSION=3 ..
  make
  cp ./Oscillator.fmu ../..
  cd ../../../mass-right-fmi
fi

python3 -m venv ../sover-fmi/.venv
. ../sover-fmi/.venv/bin/activate
pip install -r requirements.txt

fmiprecice fmi-settings.json precice-settings.json
python3 ../solver-fmi/calculate-error.py ../mass-left-fmi/fmi-settings.json ../mass-left-fmi/precice-settings.json ../mass-right-fmi/fmi-settings.json ../mass-right-fmi/precice-settings.json Mass-Right

close_log
