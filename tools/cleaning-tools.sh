#!/usr/bin/env sh

error() {
    echo "Error: $1" >&2
    exit 1
}

clean_tutorial() {
    (
        set -e -u
        cd "$1"
        echo "# Cleaning up all cases in $(pwd)..."
        rm -rfv ./precice-run/

        # Run clean.sh if it exists in the base tutorial directory
        if test -f "clean.sh"; then
            ./clean.sh
        fi

        for case in */; do
            if [ "${case}" = images/ ] || [ "${case}" = reference-results/ ]; then
                continue
            fi
            case "${case}" in solver*)
                continue
            esac
            (cd "${case}" && ./clean.sh || echo "No cleaning script in ${case} - skipping")
        done
    )
}

clean_precice_logs() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up preCICE logs in $(pwd)"
        rm -fv ./precice-*-iterations.log \
            ./precice-*-convergence.log \
            ./precice-*-watchpoint-*.log \
            ./precice-*-watchintegral-*.log \
            ./core
        rm -rfv ./precice-profiling/ profiling.json trace.json
        rm -rfv ./precice-exports/
    )
}

clean_case_logs() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up general case logs in $(pwd)"
        CASENAME="$(readlink -f "$0" | xargs dirname | xargs basename)"
        rm -fv "./$CASENAME.log"
    )
}

clean_calculix() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up CalculiX case in $(pwd)"
        rm -fv ./*.cvg ./*.dat ./*.frd ./*.sta ./*.12d ./*.rout spooles.out dummy
        rm -fv WarnNodeMissMultiStage.nam
        rm -fv ./*.eig
        rm -fv ./*.vtk
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_codeaster() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up code_aster case in $(pwd)"
        rm -fv ./*.mess ./*.resu ./*.rmed
        rm -rfv ./REPE_OUT/*
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_dealii() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up deal.II case in $(pwd)"
        rm -rfv ./dealii-output/
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_fenics() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up FEniCS case in $(pwd)"
        rm -rfv ./output/
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_nutils() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up Nutils case in $(pwd)"
        rm -fv ./*.vtk
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_openfoam() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up OpenFOAM case in $(pwd)"
        if [ -n "${WM_PROJECT:-}" ] || error "No OpenFOAM environment is active."; then
            # shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
            . "${WM_PROJECT_DIR}/bin/tools/CleanFunctions"
            cleanCase > /dev/null
            rm -rfv 0/uniform/functionObjects/functionObjectProperties history
        fi
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_su2() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up SU2 case in $(pwd)"
        rm -fv ./restart_flow_*.dat ./restart_flow_*.csv forces_breakdown.dat ./surface_flow_*.csv ./flow_*.vtk ./history_*.vtk ./history.vtk ./history_*.csv ./history.csv ./surface_flow_*.vtu ./flow_*.vtu
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_aste() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up ASTE results in $(pwd)"
        rm -fv result.vtk result.stats.json
        rm -fvr fine_mesh coarse_mesh mapped
    )
}

clean_dune() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up DUNE case in $(pwd)"
        rm -fv ./dgfparser.log
        rm -fv ./*.pvd
        rm -fv ./*.vtu
        rm -rfv ./output/
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_dumux() {
   (
        set -e -u
        cd "$1"
        echo "- Cleaning up DuMuX case in $(pwd)"
        rm -fv ./*.vtu
        rm -fv ./*.pvd
        clean_precice_logs .
        clean_case_logs .
   )
}

clean_fmi() {
    (
        set -e -u
        cd "$1"
        echo "- Cleaning up FMI case in $(pwd)"
        rm -rfv ./output/
        clean_precice_logs .
        clean_case_logs .
    )
}

clean_matlab(){
    (
	set -e -u
	cd "$1"
	echo "- Cleaning up MATLAB case in $(pwd)"
	clean_precice_logs .
	clean_case_logs .
    )	
}
