#!/bin/sh
set -e -u

python3 ./generate_mesh.py > all.msh
ccx_preCICE -i solid -precice-participant Solid
