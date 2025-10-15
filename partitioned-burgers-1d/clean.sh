#!/usr/bin/env sh
set -e -u

rm -f solver-scipy/full_domain.npz
rm -f images/full_domain_evolution.png images/full_domain_timestep_slice.png images/gradient_timestep_slice.png images/initial_condition.png
rm -f initial_condition.npz # comment out?