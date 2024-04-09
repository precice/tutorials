#!/bin/bash
# Run this script at the root of the repository to generate PNG files from each precice-config.xml

set -e -u

visualize_config(){
  echo "Visualizing the configuration in $1"
  (
    cd "$1"
    if [ -f precice-config.xml ]; then
      mkdir -p images
      precice-config-visualizer precice-config.xml | dot -Tpng > "images/tutorials-$1-precice-config.png"
    fi
  )
}

export -f visualize_config

IGNORE="partitioned-heat-conduction-direct"
tutorials=$(find . -maxdepth 1 -type d -not -name ".*" | grep -vE $IGNORE | sed "s/^.\///")

parallel visualize_config ::: "$tutorials"
