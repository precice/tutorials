#!/bin/sh
set -e -u

. ../../tools/log.sh

if [ ! -f PIDcontroller.fmu ]; then
  log cd fmu
  log rm -rf build
  log mkdir build
  log cd build
  log cmake -DFMI_TYPE=CS -DFMI_VERSION=3 ..
  log make
  log cp ./PIDcontroller.fmu ../..
  log cd ../../
fi

log fmiprecice ./fmi-settings.json ./precice-settings.json

close_log
