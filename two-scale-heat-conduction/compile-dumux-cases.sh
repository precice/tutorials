#!/usr/bin/env sh
set -e -u

cd macro-dumux/build-cmake/appl
make macro_dumux
cd ../../../micro-dumux/build-cmake/appl
make
cd ../../../

# Move macro-dumux and micro-dumux executables to the participant folder level
mv macro-dumux/build-cmake/appl/macro_dumux macro-dumux/
mv micro-dumux/build-cmake/appl/micro_sim.cpython-310-x86_64-linux-gnu.so micro-dumux/
