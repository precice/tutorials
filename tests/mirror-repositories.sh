#!/usr/bin/env bash
set -u

working_directory="${1:thirdparty}"
mkdir -p "$working_directory"
cd "$working_directory" || exit 1

repositories=(
  https://bitbucket.org/davidscn/mercurydpm.git
  https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git
  https://git.iws.uni-stuttgart.de/dumux-appl/dumux-phasefield.git
  https://gitlab.dune-project.org/core/dune-common.git
  https://gitlab.dune-project.org/core/dune-geometry.git
  https://gitlab.dune-project.org/core/dune-grid.git
  https://gitlab.dune-project.org/core/dune-istl.git
  https://gitlab.dune-project.org/core/dune-localfunctions.git
  https://gitlab.dune-project.org/extensions/dune-foamgrid.git
  https://gitlab.dune-project.org/extensions/dune-subgrid.git
  https://gitlab.dune-project.org/staging/dune-functions.git
  https://gitlab.dune-project.org/staging/dune-typetree.git
  https://gitlab.dune-project.org/staging/dune-uggrid.git
  https://gitlab.dune-project.org/extensions/dune-SPGrid.git
)

echo "Mirroring ${#repositories[@]} git repositories into $working_directory:"

counter=0
errors=0
for url in "${repositories[@]}"; do
  (( counter++ ))
  target_dir="${url#https://}"
  if [ ! -d "$target_dir" ]; then
    echo "${counter}. Cloning repository $url..."
    if ! git clone --mirror -q --progress "$url" "$target_dir"
    then
      (( errors++ ))
    fi
  else
    echo "${counter}. Updating repository $url..."
    cd "$target_dir" || exit 1
    if ! git fetch --tags
    then
      (( errors++ ))
    fi
    cd - > /dev/null || exit 1
  fi
done

if [ "${errors}" -eq 0 ]; then
  echo "Mirroring complete."
else
  echo "Mirroring incomplete (${errors} repositories failed). Try running this script again."
fi
