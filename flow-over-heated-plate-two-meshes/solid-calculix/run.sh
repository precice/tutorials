#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 ./generate_mesh.py > all.msh
log ccx_preCICE -i solid -precice-participant Solid

close_log
