#!/bin/bash
set -e -u

. ../tools/cleaning-tools.sh

clean_abaqus .

rm -rfv output/
rm -fv *.log
rm -rfv __pycache__
rm -fv .nfs*
rm -rfv ruc_*/
