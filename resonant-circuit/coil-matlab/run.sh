#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

# Run MATLAB code without GUI
matlab -nodisplay -r "coil;exit;"

close_log