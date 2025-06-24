#!/bin/sh
set -e -u

rm ./*.png ./*.mat

. ../../tools/cleaning-tools.sh

clean_matlab .
