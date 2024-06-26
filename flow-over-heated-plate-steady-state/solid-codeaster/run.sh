#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

echo "Warning: this case requires a manual preparation step for code_aster."
echo "You also need to set an absolute path as exchange-directory in precice-config.xml."
echo "See the tutorial and code_aster adapter documentation pages for more:"
echo "https://precice.org/adapter-code_aster.html"
echo ""

export TUTORIAL_ROOT="${PWD}"
export PRECICE_PARTICIPANT=Solid
as_run --run solid.export

close_log
