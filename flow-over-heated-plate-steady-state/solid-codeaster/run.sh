#!/bin/sh
set -e -u

export TUTORIAL_ROOT=${PWD}
export PRECICE_PARTICIPANT=Solid
as_run --run solid.export
