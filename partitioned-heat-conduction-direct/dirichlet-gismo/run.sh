#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

partitioned-heat-conduction -s 0 ../precice-config.xml --plot

close_log
