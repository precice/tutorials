#!/bin/sh
set -e -u

echo "Warning: this case requires a manual preparation step for code_aster."
echo "See the tutorial and code_aster adapter documentation pages for more:"
echo "https://www.precice.org/adapter-code_aster.html\n"

export TUTORIAL_ROOT=${PWD}
export PRECICE_PARTICIPANT=Solid
as_run --run solid.export
