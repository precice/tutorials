#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -f all.msh ]; then
    echo "Downloading the pre-processed Solid case..."
    wget -nv -O - https://syncandshare.lrz.de/dl/fi3c9Xt5UzBc5hJvXzsLBHXn/Solid.tar.gz | tar -xzv
    mv ./Solid/* .
    rm -r ./Solid
    sed -i 's/Solid/\./g' solid.inp
fi

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
ccx_preCICE -i solid -precice-participant Solid

close_log
