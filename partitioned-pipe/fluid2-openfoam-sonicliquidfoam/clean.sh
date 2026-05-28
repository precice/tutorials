#!/usr/bin/env sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_openfoam .

if [ -f ../precice-config.xml.orig ]; then
  echo "Restoring the precice-config.xml from precice-config.xml.orig"
  cp -r ../precice-config.xml.orig ../precice-config.xml
  rm ../precice-config.xml.orig
fi