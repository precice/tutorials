#!/bin/sh
cd ${0%/*} || exit 1                        # Run from this directory
. $WM_PROJECT_DIR/bin/tools/CleanFunctions  # Tutorial clean functions

(
    cleanAdiosOutput
    cleanAuxiliary
    cleanDynamicCode
    cleanOptimisation
    cleanPostProcessing
    cleanTimeDirectories
    cleanCase0
    rm -rf ./preCICE-output/
    rm -rf ./precice-*/
    rm -rf ../precice-*/
    rm -rf ./export/
    rm -rf ./processor*
    rm -rf ./postProcessing/*/*[1-9]*
    rm -f log.*
    rm -f *.log
    rm -rf ./*.json
)


#------------------------------------------------------------------------------
