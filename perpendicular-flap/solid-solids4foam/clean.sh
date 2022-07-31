#!/bin/bash

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cleanCase
rm -rf log.* *.log
rm -f precice-Solid-events.json
rm -rf history
