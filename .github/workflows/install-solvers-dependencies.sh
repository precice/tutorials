#!/bin/bash -e

# Install OpenFOAM v2012
wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get install openfoam2012-dev

# Install preCICE dependencies
sudo apt-get install -y cmake libeigen3-dev libxml2-dev libboost-all-dev petsc-dev python3-dev python3-numpy

# test
