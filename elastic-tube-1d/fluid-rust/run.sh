#!/bin/sh
set -e -u

. ../../tools/log.sh

log cargo run --release ../precice-config.xml 

close_log
