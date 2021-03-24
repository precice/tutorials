#!/bin/bash

# Solid participant

# Run this script in one terminal and the "runFluid" script in another terminal.
# These scripts present how the two participants would be started manually.
# Alternatively, you may execute the "Allrun" script in one terminal.

# The script "Allclean" cleans-up the result and log files.

echo "Running the Solid participant..."

# Run
ccx_preCICE -i flap -precice-participant Solid 2>&1 | tee Solid.log
