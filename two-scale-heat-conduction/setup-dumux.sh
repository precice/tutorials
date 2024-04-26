#!/bin/sh
set -e -u

# This script sets up a DUNE environment in the working directory to solve the two-scale-heat-conduction problem with DuMuX on one or both scales

# Clean any old leftover dumux or dune folders
rm -rfv dumux/ dumux-adapter/ dumux-phasefield/
rm -rfv dune-*/
rm -rfv install*

# Get the DuMuX install script and install it
wget https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/releases/3.7/bin/installdumux.py
DUNE_CONTROL_PATH=. python3 installdumux.py
# clear build directories
cd dumux
rm -r dune-common/build-cmake/dune-env/lib/dunecontrol || true
DUNE_CONTROL_PATH=. ./dune-common/bin/dunecontrol exec rm -r build-cmake
cd ..

# Take out all the module folders from the dumux/ folder and remove the dumux/ folder
mv dumux dumux-install
mv dumux-install/* ./
rm -r dumux-install

# Get additional required DUNE modules
# DuMux-preCICE adapter
git clone https://github.com/precice/dumux-adapter.git
# DuMux phasefield implementation
git clone -b cell_problems https://git.iws.uni-stuttgart.de/dumux-appl/dumux-phasefield.git
# DUNE SPGrid for periodic boundary conditions
DUNE_CONTROL_PATH=. python3 dumux/bin/installexternal.py spgrid

# Re-build environment
DUNE_CONTROL_PATH=. ./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts all

# Compile and move macro-dumux and micro-dumux executables to the participant folder level
./compile-dumux-cases.sh
