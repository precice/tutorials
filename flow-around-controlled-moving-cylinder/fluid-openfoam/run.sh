#!/bin/sh
cd ${0%/*} || exit 1                        # Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions

runApplication blockMesh
runApplication transformPoints -scale '(0.0016  0.0016 1)'
mv log.transformPoints log.transformPointsscale
runApplication transformPoints -translate '(0.0  0.0 -0.05)'
mv log.transformPoints log.transformPointsscale

restore0Dir

touch fluid-openfoam.foam

# Single
runApplication $(getApplication)

# Parallel
#runApplication decomposePar
#runParallel $(getApplication)


# ../../tools/run-openfoam.sh "$@"
# . ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs




#------------------------------------------------------------------------------
