#!/bin/bash
cd ${0%/*} || exit 1    		    # Run from this directory

# Solid participant

# Run this script in one terminal and the execute the Fluid participant 
# in another terminal.
# These scripts present how the two participants would be started manually.

for i in "$@"
do
case $i in
    -p=*|--prefix=*)
    SEARCHPATH="${i#*=}"
    shift # past argument=value
    ;;
    *)
      # unknown option
    ;;
esac
done

if test -f "elasticity"; then
    ./elasticity parameters.prm
    elif [ ! -z ${SEARCHPATH} ]; then
    ${SEARCHPATH}/elasticity parameters.prm
    else
    echo "Unable to find the executable 'elasticity'. Either specify a prefix (-p=/path/to/elasticity) or make it discoverable at runtime (e.g. export PATH)"
fi

