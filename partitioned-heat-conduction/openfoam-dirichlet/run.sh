#!/bin/sh
set -e -u

blockMesh
touch openfoam-dirichlet.foam

if test -f "heatTransfer"; then
    ./"heatTransfer"
else
	echo "Unable to find the executable 'heatTransfer'. Please compile the solver for this tutorial using the command 'wmake' in the solver directory (../openfoam-solver/) and copy the executable into this directory to run this tutorial."
fi

. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
