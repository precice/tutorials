#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

solver=./rigid_body_solver
if [ -f "${solver}" ]; then
    ${solver}
else 
    echo "Unable to locate the executable ${solver}. Have a look at the README for building instructions."
fi

close_log
