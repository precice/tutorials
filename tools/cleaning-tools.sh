#!/bin/sh
set -e -u

error() {
    echo "Error: $1" >&2; exit 1;
}

clean_tutorial() {
    (
        cd "$1"
        echo "-- Cleaning up all cases in $(pwd)..."
        rm -rfv ./precice-run/

        for case in */; do
            if [ "${case}" = images/ ]; then
                continue
            fi;
            (cd "${case}" && ./clean.sh)
        done
    )
}

clean_precice_logs() {
    (
        cd "$1"
        echo "---- Cleaning up preCICE logs in $(pwd)"
        rm -fv ./precice-*-iterations.log \
               ./precice-*-convergence.log \
               ./precice-*-events.json \
               ./precice-*-events-summary.log \
               ./precice-postProcessingInfo.log \
               ./precice-*-watchpoint-*.log \
               ./precice-*-watchintegral-*.log
    )
}

clean_calculix() {
    (
        cd "$1"
        echo "--- Cleaning up CalculiX case in $(pwd)"
        rm -fv ./*.cvg ./*.dat ./*.frd ./*.sta ./*.12d spooles.out
        clean_precice_logs .
    )
}

clean_codeaster() {
    (
        cd "$1"
        echo "--- Cleaning up code_aster case in $(pwd)"
        rm -fv ./*.mess ./*.resu ./*.rmed  # TODO: Check this list
        rm -rfv ./REPE_OUT/*
        clean_precice_logs .
    )
}

clean_dealii() {
    (
        cd "$1"
        echo "--- Cleaning up deal.II case in $(pwd)"
        rm -fv ./dealii-output/solution-*.vtk
        clean_precice_logs .
    )
}

clean_fenics() {
    (
        cd "$1"
        echo "--- Cleaning up FEniCS case in $(pwd)"
        rm -fv ./*.pvd spooles.out FSI-S/* # TODO: Check this list
        clean_precice_logs .
    )
}

clean_nutils() {
    (
        cd "$1"
        echo "--- Cleaning up Nutils case in $(pwd)"
        clean_precice_logs .
    )
}

clean_openfoam() {
    (
        cd "$1"
        echo "--- Cleaning up OpenFOAM case in $(pwd)"
        if [ -n "${WM_PROJECT:-}" ] || error "No OpenFOAM environment is active."; then
            # shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
            . "${WM_PROJECT_DIR}/bin/tools/CleanFunctions"
            cleanCase
        fi
        clean_precice_logs .
    )
}

clean_su2() {
    (
        cd "$1"
        echo "--- Cleaning up SU2 case in $(pwd)"
        rm -fv ./restart_flow_*.dat forces_breakdown.dat ./surface_flow_*.csv ./flow_*.vtk ./history_*.vtk
        clean_precice_logs .
    )
}