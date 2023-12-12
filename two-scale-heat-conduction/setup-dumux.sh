# This script sets up a DUNE environment in the working directory to solve the two-scale-heat-conduction problem with DuMuX on one or both scales

# Get the DuMuX install script and run it
wget https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/bin/installdumux.py
python3 installdumux.py

# Take out all the module folders from the dumux/ folder and remove the dumux/ folder
cd dumux
mv * ..
cd ..
rm -r dumux

# Additional required DuMuX repositories
# DuMuX-preCICE adapter
git clone -b v1.0.0 https://github.com/precice/dumux-adapter.git
cd dumux-adapter
./dune-common/bin/dunecontrol --only=dumux-precice all

cd ..

# DuMuX modules required for the micro problem
git clone https://git.iws.uni-stuttgart.de/dumux-appl/dumux-phasefield.git

python3 dumux/bin/installexternal.py spgrid

# Clear the CMake cache
./dune-common/bin/dunecontrol bexec rm -r CMakeFiles CMakeCache.txt

# Reconfigure and build DuMuX
./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts all


