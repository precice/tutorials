#!/usr/bin/env bash
set -e -u

CASENAME="$(pwd | xargs basename)"
LOGFILE="$CASENAME.log"
export LOGFILE

STARTDATE="$(date --rfc-email)"
STARTTIME="$(date +%s)"
echo "Started on: $STARTDATE" | tee "$CASENAME.log" 2>&1

close_log() {
    echo "Started on:  $STARTDATE" | tee --append "$LOGFILE" 2>&1
    ENDDATE="$(date --rfc-email)"
    ENDTIME="$(date +%s)"
    echo "Finished on: $ENDDATE" | tee --append "$LOGFILE" 2>&1
    echo "Duration:    $((ENDTIME-STARTTIME)) seconds (wall-clock time, including time waiting for participants)" | tee --append "$LOGFILE" 2>&1
}
