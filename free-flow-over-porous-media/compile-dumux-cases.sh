#!/usr/bin/env sh
set -e -u

cd free-flow-dumux/build-cmake/appl
make free_flow_dumux
cd ../../../porous-media-dumux/build-cmake/appl
make porous_media_dumux
cd ../../../

# Move free-flow-dumux and porous-media-dumux executables to the participant folder level
mv free-flow-dumux/build-cmake/appl/free_flow_dumux free-flow-dumux/
mv porous-media-dumux/build-cmake/appl/porous_media_dumux porous-media-dumux/
