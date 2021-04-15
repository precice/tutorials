#!/bin/sh
set -e -u

solver=./rigid_body_solver
if [ -f "${solver}" ]; then
    ${solver}
else 
    echo "Unable to locate the executable ${solver}. Have a look at the README for building instructions."
fi
