#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

usage() { echo "Usage: cmd [-s] [-p n] [-l <path-to-DUNE-common>]" 1>&2; exit 1; }

if [ ! -d "build-cmake" ]; then
  echo "Solver not built. Building now..."
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
    ../dune-common/bin/dunecontrol --opts=../dumux/cmake.opts --only=micro_sim all
  else
    DUNE_COMMON_PATH=$DUNE_COMMON_PATH_ARG
    export DUNE_CONTROL_PATH=$DUNE_COMMON_PATH:$CASE_DIR
    "$DUNE_COMMON_PATH"/dune-common/bin/dunecontrol --opts="$DUNE_COMMON_PATH"/dumux/cmake.opts --only=micro_sim all
  fi
else
  echo "build-cmake folder found."
  cd build-cmake
  make micro_sim
  cd ..
fi
# Move micro_sim executable to the participant folder level
mv ./build-cmake/appl/micro_sim*.so .

# Check if no input argument was provided
if [ -z "$*" ] ; then
  echo "No input argument provided. Micro Manager is launched in serial"
  micro-manager-precice micro-manager-config.json
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    micro-manager-precice micro-manager-config.json
    ;;
  p)
    mpiexec -n "$2" micro-manager-precice micro-manager-config.json
    ;;
  *)
    usage
    ;;
  esac
done

close_log
