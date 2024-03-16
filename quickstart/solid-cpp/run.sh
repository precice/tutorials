#!/bin/sh
set -e -u

. ../../tools/log.sh

solver=./rigid_body_solver
if [ -f "${solver}" ]; then
    log ${solver}
else 
    log echo "Unable to locate the executable ${solver}. Have a look at the README for building instructions."
fi

close_log
