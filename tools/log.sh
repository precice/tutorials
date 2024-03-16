#!/bin/sh
set -e -u

CASENAME="$(pwd | xargs basename)"
export CASENAME

STARTDATE="$(date --rfc-email)"
STARTTIME="$(date +%s)"
echo "Started on: $STARTDATE" | tee "$CASENAME.log" 2>&1

log() {
    "$@" | tee --append "$CASENAME.log" 2>&1
}

close_log() {
    echo "Started on:  $STARTDATE" | tee --append "$CASENAME.log" 2>&1
    ENDDATE="$(date --rfc-email)"
    ENDTIME="$(date +%s)"
    echo "Finished on: $ENDDATE" | tee --append "$CASENAME.log" 2>&1
    echo "Duration:    $((ENDTIME-STARTTIME)) seconds (wall-clock time, including time waiting for participants)" | tee --append "$CASENAME.log" 2>&1
}
