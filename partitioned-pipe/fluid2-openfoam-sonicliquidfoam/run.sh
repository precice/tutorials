#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

# The precice-config.xml has the PressureGradient field commented out.
# That improves the results for the compressible case, so we un-comment it.
if [ ! -f ../precice-config.xml.orig ]; then
  echo "Modifying the ../precice-config.xml to enable PressureGradient (see precice-config.xml.orig for the original)"
  cp -r ../precice-config.xml ../precice-config.xml.orig
  sed -i "s,<!--,,g" ../precice-config.xml
  sed -i "s,-->,,g" ../precice-config.xml
fi

blockMesh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
