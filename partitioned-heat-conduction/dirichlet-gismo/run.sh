#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

./gismo-executable -c ../precice-config.xml --plot -s 0

close_log
