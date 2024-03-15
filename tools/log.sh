#!/bin/sh
set -e -u

CASENAME="$(readlink -f "$0" | xargs dirname | xargs basename)"
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
    echo "Duration:    $((ENDTIME-STARTTIME)) seconds" | tee --append "$CASENAME.log" 2>&1
}
