#!/bin/sh
set -e -u

export OMP_NUM_THREADS=1 NUTILS_NPROCS=4

python3 transport.py
