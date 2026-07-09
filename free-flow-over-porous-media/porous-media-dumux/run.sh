#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1


usage() { echo "Usage: cmd [-l <path-to-DUNE-common>]" 1>&2; exit 1; }

DUNE_COMMON_PATH_SET=
DUNE_COMMON_PATH_ARG=

while getopts ":l:" opt; do
  case ${opt} in
  l)
    DUNE_COMMON_PATH_SET=1
    DUNE_COMMON_PATH_ARG="$OPTARG"
    ;;
  *)
    usage
    ;;
  esac
done

if [ ! -d "build-cmake" ]; then
  echo "Solver not built. Building now..."
  CASE_DIR=$(pwd)/..

  if [ -z "$DUNE_COMMON_PATH_SET" ]; then
    ../dune-common/bin/dunecontrol --opts=../dumux/cmake.opts --only=porous_media_dumux all
  else
    export DUNE_CONTROL_PATH=$DUNE_COMMON_PATH_ARG:$CASE_DIR
    "$DUNE_COMMON_PATH_ARG"/dune-common/bin/dunecontrol --opts="$DUNE_COMMON_PATH_ARG"/dumux/cmake.opts --only=porous_media_dumux all
  fi
else
  echo "build-cmake folder found."
  cd build-cmake
  make porous_media_dumux
  cd ..
fi
# Move porous_media_dumux executable to the participant folder level
mv build-cmake/solver-dumux/porous_media_dumux .

./porous_media_dumux params.input

close_log