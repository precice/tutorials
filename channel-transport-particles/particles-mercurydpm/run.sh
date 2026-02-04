#!/usr/bin/env bash
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
    "${EXE}"
    exit 0
fi

EXE="ChannelTransport"
# If it is in the global path
if [ -n "$(command -v "${EXE}")" ]; then
    "${EXE}"
    exit 0
fi

# If it has been copied to the local directory
if test -f "ChannelTransport"; then
    ./"${EXE}"
else
    echo "Unable to find the executable ${EXE}. Either specify the executable explicitly (-e=/path/to/ChannelTransport) or make it discoverable at runtime (e.g. export PATH). In the MercuryDPM build, the executable is compiled in the directory /path/to/mercurydpm/build/Drivers/PreCICE/ChannelTransport"
fi
