#!/bin/sh
set -e -u

# Remove the old directory and copy the uninitialized field
rm -rf ./0
cp -r ./0.orig 0
# Initialize the new field
funkySetFields -keepPatches -field T -expression '1+pow(pos().x,2)+(3*pow(pos().y,2))+1.2*time()' -time '0'
