#!/bin/bash
# Run this script at the root of the repository to generate PNG files from each precice-config.xml

set -e -u

tutorials=$(find . -maxdepth 1 -type d -not -name ".*" | sed "s/^.\///")

for tutorial in $tutorials; do
  if [ "${tutorial}" = tools ] || [ "${tutorial}" = changelog-entries ] || [ "${tutorial}" = aste-turbine ] || [ "${tutorial}" = partitioned-heat-conduction-direct ] || [ "${tutorial}" = flow-over-heated-plate-nearest-projection ] || [ "${tutorial}" = flow-over-heated-plate-two-meshes ] || [ "${tutorial}" = partitioned-heat-conduction-complex ]; then
    continue
  else
    echo "Visualizing the configuration in $tutorial"
    (
      cd "$tutorial"
      precice-config-visualizer precice-config.xml | dot -Tpng > "images/tutorials-$tutorial-precice-config.png"
    )
  fi
done
