#!/usr/bin/env sh
set -e -u

# To compile macro-dumux and micro-dumux from scratch or recompile them after changes
./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts --only=macro_dumux all
./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts --only=micro_sim all

# Alternatively, you can manually recompile macro-dumux and micro-dumux when the `build-cmake` folder is present by uncommenting the following lines:
# (
# cd macro-dumux/build-cmake/appl
# make macro_dumux
# )
# (
# cd micro-dumux/build-cmake/appl
# make micro_sim
# )

# Move macro-dumux and micro-dumux executables to the participant folder level
mv macro-dumux/build-cmake/appl/macro_dumux macro-dumux/
mv micro-dumux/build-cmake/appl/micro_sim.cpython-*.so micro-dumux/
