#!/bin/sh
set -e -u

CASENAME="$(readlink -f "$0" | xargs dirname | xargs basename)"
export CASENAME
