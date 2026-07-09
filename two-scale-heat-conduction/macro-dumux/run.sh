#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

usage() { echo "Usage: cmd [-s] [-p n] [-l <path-to-DUNE-common>]" 1>&2; exit 1; }

SERIAL_RUN=
PARALLEL_RUN=
RANK_COUNT=
DUNE_COMMON_PATH_SET=
DUNE_COMMON_PATH_ARG=

while getopts ":sp:l:" opt; do
  case ${opt} in
  s)
    SERIAL_RUN=1
    ;;
  p)
    RANK_COUNT="$OPTARG"
    PARALLEL_RUN=1
    ;;
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
    ../dune-common/bin/dunecontrol --opts=../dumux/cmake.opts --only=macro_dumux all
  else
    export DUNE_CONTROL_PATH=$DUNE_COMMON_PATH_ARG:$CASE_DIR
    "$DUNE_COMMON_PATH_ARG"/dune-common/bin/dunecontrol --opts="$DUNE_COMMON_PATH_ARG"/dumux/cmake.opts --only=macro_dumux all
  fi
else
  echo "build-cmake folder found."
  cd build-cmake
  make macro_dumux
  cd ..
fi
# Move macro-dumux executable to the participant folder level
mv ./build-cmake/appl/macro_dumux .

if [ -n "$SERIAL_RUN" ] && [ -n "$PARALLEL_RUN" ]; then
  echo "Cannot run both serial and parallel. Choose one option."
  usage
elif [ -z "$SERIAL_RUN" ] && [ -z "$PARALLEL_RUN" ]; then
  echo "No run option provided. The macro solver is launched in serial."
fi

if [ -n "$PARALLEL_RUN" ]; then
  mpiexec -n "$RANK_COUNT" macro_dumux params.input
else
  ./macro_dumux params.input
fi

close_log
