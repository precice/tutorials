#!/bin/bash

echo "Cleaning..."

. $WM_PROJECT_DIR/bin/tools/CleanFunctions

# Participant 1: Solid
Participant1="solid-codeaster"
cd ${Participant1}
    rm -fv ${Participant1}/solid.mess
    rm -fv ${Participant1}/solid.resu
    rm -fv ${Participant1}/solid.rmed
    rm -fv ${Participant1}.log
    rm -fvr ${Participant1}/REPE_OUT
    mkdir ${Participant1}/REPE_OUT
    rm -fv \
      precice-*.log \
      precice-*-events.json
cd ..

# Participant 2: Fluid
Participant2="Fluid"
rm -fv ${Participant2}.log
cd ${Participant2}
  cleanCase
  touch Fluid.foam
  # Remove the log files
  rm -fv ${Participant2}_blockMesh.log
  rm -fv ${Participant2}_checkMesh.log
  rm -fv ${Participant2}_decomposePar.log
  rm -fv ${Participant2}_reconstructPar.log
  rm -fv \
    precice-*.log \
    precice-*-events.json
cd ..



# Output files for preCICE versions before 1.2:
rm -fv \
    iterations-${Participant1}.txt iterations-${Participant2}.txt \
    convergence-${Participant1}.txt convergence-${Participant2}.txt \
    Events-${Participant1}.log Events-${Participant2}.log \
    EventTimings-${Participant1}.log EventTimings-${Participant2}.log

# Remove the preCICE address files
rm -rfv precice-run
rm -fv .*.address

echo "Cleaning complete!"
