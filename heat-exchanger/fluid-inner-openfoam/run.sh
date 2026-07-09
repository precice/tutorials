#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -f constant/polyMesh/boundary ]; then
  echo "Downloading and extracting the Inner-Fluid mesh..."
  wget -nv -O - https://syncandshare.lrz.de/dl/fiNsYGC1DKzgio4jS5NhsXg7/polyMesh.org.tar.gz | tar -xzv -C constant
  mv constant/polyMesh.org constant/polyMesh
fi

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
