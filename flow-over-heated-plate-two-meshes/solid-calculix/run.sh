#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

python3 ./generate_mesh.py > all.msh
ccx_preCICE -i solid -precice-participant Solid

close_log
