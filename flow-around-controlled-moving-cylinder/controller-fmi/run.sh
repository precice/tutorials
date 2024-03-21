#!/bin/bash
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

fmiprecice ./fmi-settings.json ./precice-settings.json

close_log
