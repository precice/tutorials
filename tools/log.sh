#!/usr/bin/env bash
set -e -u

CASENAME="$(pwd | xargs basename)"
LOGFILE="$CASENAME.log"
export LOGFILE

rfc_email_date() {
    date --rfc-email 2>/dev/null || date '+%a, %d %b %Y %H:%M:%S %z'
}

STARTDATE="$(rfc_email_date)"
STARTTIME="$(date +%s)"
echo "Started on: $STARTDATE" | tee "$CASENAME.log" 2>&1

close_log() {
    echo "Started on: $STARTDATE" | tee -a "$LOGFILE" 2>&1
    ENDDATE="$(rfc_email_date)"
    ENDTIME="$(date +%s)"
    echo "Finished on: $ENDDATE" | tee -a "$LOGFILE" 2>&1
    echo "Duration:    $((ENDTIME-STARTTIME)) seconds (wall-clock time, including time waiting for participants)" | tee -a "$LOGFILE" 2>&1
}
