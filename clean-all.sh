#!/bin/sh
set -e -u

echo "- Cleaning up all tutorials..."

find . -maxdepth 2 -mindepth 2 -name clean-tutorial.sh -execdir sh -c './clean-tutorial.sh' \;