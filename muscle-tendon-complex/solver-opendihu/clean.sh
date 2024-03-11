#!/bin/bash
set -e -u

rm -rf ./*.log
rm -rf build_release
rm -rf precice-profiling

# remove executables from partipant folders
rm -f ../muscle-opendihu/muscle-solver
rm -f ../tendon-bottom-opendihu/tendon-solver
rm -f ../tendon-top-A-opendihu/tendon-solver
rm -f ../tendon-top-B-opendihu/tendon-solver
