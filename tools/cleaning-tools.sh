#!/bin/sh

error() {
    echo "Error: $1" >&2
    exit 1
}

clean_tutorial() {
    (
        set -e -u
        cd "$1"
        echo "-- Cleaning up all cases in $(pwd)..."
        rm -rfv ./precice-run/

        for case in */; do
            if [ "${case}" = images/ ]; then
                continue
            fi
            (cd "${case}" && ./clean.sh || echo "No cleaning script in ${case} - skipping")
        done
    )
}

clean_precice_logs() {
    (
        set -e -u
        cd "$1"
        echo "---- Cleaning up preCICE logs in $(pwd)"
        rm -fv ./precice-*-iterations.log \
            ./precice-*-convergence.log \
            ./precice-*-events.json \
            ./precice-*-events-summary.log \
            ./precice-postProcessingInfo.log \
            ./precice-*-watchpoint-*.log \
            ./precice-*-watchintegral-*.log \
            ./core
    )
}

clean_calculix() {
    (
        set -e -u
        cd "$1"
        echo "--- Cleaning up CalculiX case in $(pwd)"
        rm -fv ./*.cvg ./*.dat ./*.frd ./*.sta ./*.12d spooles.out dummy
        clean_precice_logs .
    )
}

clean_codeaster() {
    (
        set -e -u
        cd "$1"
        echo "--- Cleaning up code_aster case in $(pwd)"
        rm -fv ./*.mess ./*.resu ./*.rmed
        rm -rfv ./REPE_OUT/*
        clean_precice_logs .
    )
}

clean_dealii() {
    (
        set -e -u
        cd "$1"
        echo "--- Cleaning up deal.II case in $(pwd)"
        rm -rfv ./dealii-output/
        clean_precice_logs .
    )
}

clean_fenics() {
    (
        set -e -u
        cd "$1"
        echo "--- Cleaning up FEniCS case in $(pwd)"
        rm -fv ./*.pvd spooles.out FSI-S/*
        rm -rfv ./out/
        rm -rfv ./preCICE-output/
        clean_precice_logs .
    )
}

clean_nutils() {
    (
        set -e -u
        cd "$1"
        echo "--- Cleaning up Nutils case in $(pwd)"
        rm -fv ./*.vtk
        rm -rfv ./preCICE-output/
        clean_precice_logs .
    )
}

clean_openfoam() {
    (
        set -e -u
        cd "$1"
        echo "--- Cleaning up OpenFOAM case in $(pwd)"
        if [ -n "${WM_PROJECT:-}" ] || error "No OpenFOAM environment is active."; then
            # shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
            . "${WM_PROJECT_DIR}/bin/tools/CleanFunctions"
            cleanCase
            rm -rfv 0/uniform/functionObjects/functionObjectProperties
        fi
        rm -rfv ./preCICE-output/
        clean_precice_logs .
    )
}

clean_su2() {
    (
        set -e -u
        cd "$1"
        echo "--- Cleaning up SU2 case in $(pwd)"
        rm -fv ./restart_flow_*.dat forces_breakdown.dat ./surface_flow_*.csv ./flow_*.vtk ./history_*.vtk
        clean_precice_logs .
    )
}
