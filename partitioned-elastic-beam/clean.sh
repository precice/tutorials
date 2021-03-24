#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

echo "Cleaning..."


Participant1="dirichlet-calculix"
cd ${Participant1}
    # Clean the case
    rm -fv *.log
    rm -fv beam1.cvg
    rm -fv beam1.dat
    rm -fv beam1.frd
    rm -fv beam1.sta
    rm -fv beam1.12d
    rm -fv spooles.out
    rm -fv precice-*.log
    rm -fv precice-*-events.json
cd ..

Participant2="neumann-calculix"
cd ${Participant2}
    # Clean the case
    rm -fv *.log
    rm -fv beam2.cvg
    rm -fv beam2.dat
    rm -fv beam2.frd
    rm -fv beam2.sta
    rm -fv beam2.12d
    rm -fv spooles.out
    rm -fv precice-*.log
    rm -fv precice-*-events.json
cd ..
# Remove the preCICE address files
rm -rfv precice-run

echo "Cleaning complete!"
#------------------------------------------------------------------------------
