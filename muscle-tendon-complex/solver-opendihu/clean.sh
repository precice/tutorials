#!/bin/bash
set -e -u

rm -r .* ./*.log
rm -r build_release
rm -r precice-profiling

# remove executables from partipant folders
rm ../muscle-opendihu/muscle-solver
rm ../tendon-bottom-opendihu/tendon-solver
rm ../tendon-top-A-opendihu/tendon-solver
rm ../tendon-top-B-opendihu/tendon-solver
