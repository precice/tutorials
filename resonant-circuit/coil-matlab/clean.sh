#!/bin/sh
set -e -u

find . -not -name "clean.sh" -not -name "run.sh" -not -name "Solver_Coil.m" -delete
