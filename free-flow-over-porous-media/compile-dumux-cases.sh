#!/usr/bin/env sh
set -e -u

# To compile free-flow-dumux and porous-media-dumux from scratch or recompile them after changes
./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts --only=free_flow_dumux all
./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts --only=porous_media_dumux all

# Alternatively, you can manually recompile free-flow-dumux and porous-media-dumux when the `build-cmake` folder is present by uncommenting the following lines:
# (
# cd free-flow-dumux/build-cmake/solver-dumux
# make free_flow_dumux
# )

# (
# cd porous-media-dumux/build-cmake/solver-dumux
# make porous_media_dumux
# )

# Move free-flow-dumux and porous-media-dumux executables to the participant folder level
mv free-flow-dumux/build-cmake/solver-dumux/free_flow_dumux free-flow-dumux/
mv porous-media-dumux/build-cmake/solver-dumux/porous_media_dumux porous-media-dumux/
