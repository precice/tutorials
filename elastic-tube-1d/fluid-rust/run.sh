#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

pkg-config --libs --cflags libprecice
cargo run --release ../precice-config.xml 

close_log
