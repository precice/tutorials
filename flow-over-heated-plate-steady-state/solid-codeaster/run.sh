#!/bin/sh
set -e -u

. ../../tools/log.sh

log echo "Warning: this case requires a manual preparation step for code_aster."
log echo "You also need to set an absolute path as exchange-directory in precice-config.xml."
log echo "See the tutorial and code_aster adapter documentation pages for more:"
log echo "https://precice.org/adapter-code_aster.html"
log echo ""

export TUTORIAL_ROOT="${PWD}"
export PRECICE_PARTICIPANT=Solid
log as_run --run solid.export

close_log
