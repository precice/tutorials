#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

echo "Free flow solver is launched in serial."
./free_flow_dumux params.input

close_log