#!/bin/sh
set -e -u

. ../../tools/log.sh

if [ ! -f PIDcontroller.fmu ]; then
  cd fmu
  rm -rf build
  mkdir build
  cd build
  log cmake -DFMI_TYPE=CS -DFMI_VERSION=3 ..
  log make
  cp ./PIDcontroller.fmu ../..
  cd ../../
fi

log fmiprecice ./fmi-settings.json ./precice-settings.json

close_log
