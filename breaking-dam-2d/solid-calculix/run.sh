#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

ccx_preCICE -i flap -precice-participant Solid

close_log
