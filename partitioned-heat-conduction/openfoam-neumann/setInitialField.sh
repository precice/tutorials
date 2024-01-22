#!/bin/sh
set -e -u

# Remove the old directory and copy the uninitialized field
rm -rf ./0
cp -r ./0.orig 0

