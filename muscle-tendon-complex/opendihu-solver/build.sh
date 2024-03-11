#!/bin/bash

if [ -n "$OPENDIHU_HOME" ]
then 
    "$OPENDIHU_HOME/scripts/shortcuts/sr.sh" && "$OPENDIHU_HOME/scripts/shortcuts/mkorn.sh"
    # copy executables to partipant folders
    cp build_release/muscle-solver ../muscle-opendihu/
    cp build_release/tendon-solver ../tendon-bottom-opendihu/
    cp build_release/tendon-solver ../tendon-top-A-opendihu/
    cp build_release/tendon-solver ../tendon-top-B-opendihu/
else
    echo "OPENDIHU_HOME is not defined"
fi
