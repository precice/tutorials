#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

perpendicular-flap-vertex-gismo -c ../precice-config.xml --plot -r 2 

close_log
