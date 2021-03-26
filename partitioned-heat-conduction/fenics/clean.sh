#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

rm -rf out 
rm -f *.log
rm -f *-events.json
