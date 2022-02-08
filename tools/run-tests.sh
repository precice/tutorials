#!/bin/bash
set -u -e

# Execute this script from the parent directory, as ./tools/run-tests.sh

export TAG_OPENFOAM_ADAPTER=latest

cd flow-over-heated-plate/tests && MY_UID="$(id -u)" MY_GID="$(id -g)" docker-compose up
