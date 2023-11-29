#!/bin/sh
set -e -u

. ../../../tools/cleaning-tools.sh

clean_abaqus .

rm -rfv desaii*
rm -fv ./*_ruc.log
rm -fv ./*_nm1.inp
rm -fv stresses.txt
rm -fv log_ruc_run.log
