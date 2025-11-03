#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

echo "Porous medium solver is launched in serial."
./porous_media_dumux params.input

close_log