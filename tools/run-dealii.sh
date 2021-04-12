#!/bin/sh
set -e -u

EXE=""

for i in "$@"; do
    case $i in
    -e=* | --exec=*)
        EXE="${i#*=}"
        shift # past argument=value
        ;;
    *)
        # unknown option
        ;;
    esac
done

# If the executable has been defined
if [ -n "${EXE}" ]; then
    "${EXE}" parameters.prm
    exit 0
fi

EXE="elasticity"
# If it is in the global path
if [ -n "$(command -v "${EXE}")" ]; then
    "${EXE}" parameters.prm
    exit 0
fi

# If it has been copied to the local directory
if test -f "elasticity"; then
    ./"${EXE}" parameters.prm
else
    echo "Unable to find the executable ${EXE}. Either specify the executable explicitly (-e=/path/to/elasticity) or make it discoverable at runtime (e.g. export PATH)"
fi
