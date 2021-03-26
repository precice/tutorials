#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

./clean_Dirichlet.sh
./clean_Neumann.sh
rm -f *.log
rm -f *-events.json
