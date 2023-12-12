# This script sets up a DUNE environment in the subdirectory "dumux" to solve the two-scale-heat-conduction problem with DuMuX on one or both scales

# Get the DuMuX install scripts
wget https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/bin/installdumux.py
wget https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/bin/installexternal.py

# Get additional required DUNE modules
cd dumux
# DuMux-preCICE adapter
git clone -b develop https://github.com/precice/dumux-adapter.git
# DuMux phasefield implementation
git clone -b cell_problems https://git.iws.uni-stuttgart.de/dumux-appl/dumux-phasefield.git
# DUNE SPGrid for periodic boundary conditions
mkdir dune-common
python3 ../installexternal.py spgrid
rm -r dune-common
cd ..

# Install and build DuMux
python3 installdumux.py
