#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

usage() { echo "Usage: cmd [-s] [-p n] [-l <path-to-DUNE-common>]" 1>&2; exit 1; }

if [ ! -d "build-cmake" ]; then
  echo "No build-cmake folder found. Compiling macro-dumux."
  CASE_DIR=$(pwd)/..

  while getopts ":l:" opt; do
    case ${opt} in
    l)
      DUNE_COMMON_PATH_ARG=$OPTARG
      ;;
    *)
      usage
      ;;
    esac
  done
  if [ -z "$DUNE_COMMON_PATH_ARG" ]; then
    ../dune-common/bin/dunecontrol --opts=../dumux/cmake.opts --only=macro_dumux all
  else
    DUNE_COMMON_PATH=$DUNE_COMMON_PATH_ARG
    export DUNE_CONTROL_PATH=$DUNE_COMMON_PATH:$CASE_DIR
    "$DUNE_COMMON_PATH"/dune-common/bin/dunecontrol --opts="$DUNE_COMMON_PATH"/dumux/cmake.opts --only=macro_dumux all
  fi
else
  echo "build-cmake folder found."
  cd build-cmake
  make macro_dumux
  cd ..
fi
# Move macro-dumux executable to the participant folder level
mv ./build-cmake/appl/macro_dumux .

# Check if no input argument was provided
if [ -z "$*" ] ; then
  echo "No input argument provided. Macro solver is launched in serial"
  ./macro_dumux params.input
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    ./macro_dumux params.input
    ;;
  p)
    mpiexec -n "$2" macro_dumux params.input
    ;;
  *)
    usage
    ;;
  esac
done

close_log
