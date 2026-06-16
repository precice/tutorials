#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

blockMesh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

# System tests: Keep the container and the respective network alive till the end.
if [[ -v PRECICE_TUTORIALS_TESTING ]]; then
    echo "Waiting for the Particles participant to finish..."
    if [ ! -f "../particles-participant-finished.log" ]; then
        inotifywait -e create,modify,attrib --include '/particles-participant-finished\.log$' -qq ..
    fi
fi

close_log
