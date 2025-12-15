#!/usr/bin/env sh
set -e -u

# shellcheck disable=SC1091
. ../tools/cleaning-tools.sh

clean_tutorial .
clean_precice_logs .
rm -fv ./*.log
rm -fv ./*.vtu

rm -f solver-scipy/full_domain.npz
rm -f dirichlet-scipy/dirichlet.npz
rm -f neumann-scipy/neumann.npz
rm -rf output/
rm -f initial_condition.npz
