#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -f constant/polyMesh/boundary ]; then
  echo "Downloading and extracting the Outer-Fluid mesh..."
  wget -nv -O - https://syncandshare.lrz.de/dl/fiEZRQ8rcVWRkoyZvANim1R1/polyMesh.org.tar.gz | tar -xzv -C constant
  mv constant/polyMesh.org constant/polyMesh
fi

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
