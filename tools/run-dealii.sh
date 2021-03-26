#!/bin/sh
set -e -u
# Solid participant

# Run this script in one terminal and the execute the Fluid participant 
# in another terminal.
# These scripts present how the two participants would be started manually.

SEARCHPATH=""
SOLVER="elasticity"

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

# If a prefix has been defined
if [ -n "${SEARCHPATH}"  ]; then
    "${SEARCHPATH}"/"${SOLVER}" parameters.prm
fi

# If it is in the global path
if [ -n "$(command -v "${SOLVER}")" ]; then
   "${SOLVER}" parameters.prm
fi

# If it has been copied to the local directory
if test -f "elasticity"; then
    ./"${SOLVER}" parameters.prm
    else
    echo "Unable to find the executable ${SOLVER}. Either specify a prefix (-p=/path/to/elasticity) or make it discoverable at runtime (e.g. export PATH)"
fi

