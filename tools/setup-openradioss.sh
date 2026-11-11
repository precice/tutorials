#!/usr/bin/env bash

# Script that prepares an environment with aliases used in the OpenRadioss preCICE tutorials.
# This script was originally prepared for OpenRadioss, but also works with OpenCourant (as of November 2026).

# Set OPENRADIOSS_PATH in your ~/.bashrc, or edit the next line
# export OPENRADIOSS_PATH="/home/user/OpenRadioss"

# As recommended by the installation guide
export RAD_CFG_PATH=$OPENRADIOSS_PATH/hm_cfg_files
export LD_LIBRARY_PATH="$OPENRADIOSS_PATH/extlib/hm_reader/linux64/:$OPENRADIOSS_PATH/extlib/h3d/lib/linux64/:$LD_LIBRARY_PATH"
export OMP_STACKSIZE=400m

export PATH="$OPENRADIOSS_PATH/starter/cbuild_starter_linux64_gf:$OPENRADIOSS_PATH/engine/cbuild_engine_linux64_gf_ompi_precice/:$PATH"
alias radioss_starter='starter_linux64_gf'
alias radioss_engine='engine_linux64_gf_ompi_precice'

# These tools are not included in the basic installation.
# For OpenCourant, find them in https://github.com/OpenCourant/Tools/tree/main/output_converters
alias radiossVTK='$OPENRADIOSS_PATH/exec/anim_to_vtk_linux64_gf'
alias radiossCSV='$OPENRADIOSS_PATH/exec/th_to_csv_linux64_gf'

# Run `toVTK modelA` to convert all animation files modelAA* to VTK format
toVTK() {

    if [ $# -lt 1 ]; then
        echo "ERROR - Usage: toVTK ROOTNAME"
    fi

    Rootname="$1"

    for file in "${Rootname}"A*; do
        # skip if no matching files
        [ -e "$file" ] || { echo "ERROR: No files matching ${Rootname}A*"; }

        animation_number="${file#"${Rootname}A"}"

        radiossVTK "$file" > "${Rootname}_${animation_number}.vtk"

        echo "Converted ${Rootname}A${animation_number} → ${Rootname}_${animation_number}.vtk"
    done
}

# Run `toCSV modelA` to convert all animation files modelAT* to CSV format
toCSV() {

    if [ $# -lt 1 ]; then
        echo "ERROR - Usage: toVTK ROOTNAME"
    fi

    Rootname="$1"

    for file in "${Rootname}"T*; do
        # skip if no matching files
        [ -e "$file" ] || { echo "ERROR: No files matching ${Rootname}T*"; }

        history_number="${file#"${Rootname}T"}"

        radiossCSV "$file" > "${Rootname}_${history_number}.csv"

        echo "Converted ${Rootname}T${history_number} → ${Rootname}_${history_number}.csv"
    done
}