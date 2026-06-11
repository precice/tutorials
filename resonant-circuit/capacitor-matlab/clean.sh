#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_matlab .
rm -f ./*.png ./*.mat
