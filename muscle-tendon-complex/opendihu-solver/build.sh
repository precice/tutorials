#!/bin/bash

if [ -n $OPENDIHU_HOME ]
then 
    alias sr='$OPENDIHU_HOME/scripts/shortcuts/sr.sh'
    alias mkorn='$OPENDIHU_HOME/scripts/shortcuts/mkorn.sh'
    mkorn && sr
    # copy executables to partipant folders
    cp build_release/muscle-solver ../muscle-opendihu/
    cp build_release/tendon-solver ../tendon-bottom-opendihu/
    cp build_release/tendon-solver ../tendon-top-A-opendihu/
    cp build_release/tendon-solver ../tendon-top-B-opendihu/
else
    echo "OPENDIHU_HOME is not defined"
fi