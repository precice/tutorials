#!/bin/sh
set -e -u

ccx_preCICE frequency
mv frequency.eig flap.eig
ccx_preCICE -i flap -precice-participant Solid
