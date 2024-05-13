#!/usr/bin/env bash
# Run this script at the root of the repository to generate PNG files from each precice-config.xml

set -e -u

visualize_config(){
  (
    if [[ "$1" == quickstart ]]; then
      outfile="images/quickstart-precice-config"
    else
      outfile="images/tutorials-$1-precice-config"
    fi

    cd "$1"
    if [ -f precice-config.xml ]; then
      echo "Visualizing the configuration in $1"
      mkdir -p images
      precice-config-visualizer -o "$outfile.dot" precice-config.xml
      
      # Special case, to be removed once bug https://github.com/precice/config-visualizer/issues/22
      if [[ "$1" == partitioned-heat-conduction-direct ]]; then
         sed 's/compound=True;//' --in-place "$outfile.dot"
      fi
      
      dot -Tpng "$outfile.dot" > "$outfile.png"
      rm "$outfile.dot"
    fi
  )
}

export -f visualize_config

python3 -m venv .venv
. .venv/bin/activate
pip install precice-config-visualizer

tutorials=$(find . -maxdepth 1 -type d -not -name ".*" | sed "s/^.\///")

if command -v parallel &> /dev/null; then
  parallel visualize_config ::: "$tutorials"
else
  visualize_config ::: "$tutorials"
fi
