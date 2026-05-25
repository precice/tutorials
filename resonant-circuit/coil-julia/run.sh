#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

julia --project=Project.toml -e 'using Pkg; Pkg.add("DifferentialEquations"); Pkg.add("PreCICE"); Pkg.instantiate();'
julia --project=Project.toml coil.jl

close_log