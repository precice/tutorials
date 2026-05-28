#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

# The precice-config.xml has the PressureGradient field commented out.
# That improves the results for the compressible case, so we un-comment it.
cp -r ../precice-config.xml ../precice-config.xml.orig
sed -i "s,\<!--,,g" ../precice-config.xml
sed -i "s,\-->,,g" ../precice-config.xml

blockMesh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
