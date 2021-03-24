#!/bin/sh
set -e -u

echo "- Cleaning up all tutorials..."

for tutorial in */; do
    if [ "${tutorial}" = tools/ ] || [ "${tutorial}" = elastic-tube-1d/ ] || [ "${tutorial}" = partitioned-heat-conduction/ ] || [ "${tutorial}" = partitioned-heat-conduction-complex/ ] || [ "${tutorial}" = quickstart/ ]; then
    # TODO: Update when the excluded cases are also ported to the new structure
        continue
    fi;
    (cd "${tutorial}" && ./clean-tutorial.sh)
done
