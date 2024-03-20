#!/bin/sh
set -e -u

. ../../tools/log.sh

log python3 solid.py ../precice-config.xml

close_log
