#!/bin/sh

# This script cannot be used as-is and is meant to generate the picture
# images/tutorials-perpendicular-flap-displacement-all-watchpoints.png
# rendered on https://precice.org/tutorials-perpendicular-flap.html.
#
# It plots watchpoints of specific combinations of Fluid and Solid participants,
# stored in `images/`. To use this script:
# 1. For each combination you want to plot: 
#    1. Run the tutorial with that combination
#    2. Copy the respective watchpoint file to watchpoints/, with names <fluid>-<solid>.log
#    3. Clean the tutorial
# 2. Edit the script to plot the files you want, with the corresponding titles,
#    adding one line per combination.
# 3. Call the script with ./plot-all-displacements.sh, from this directory.

gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'x-displacement of the flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'x-displacement [m]'
	set term pngcairo enhanced size 900,654
	set output "images/tutorials-perpendicular-flap-displacement-all-watchpoints.png"
	plot "watchpoints/openfoam-calculix.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
	     "watchpoints/openfoam-dealii.log" using 1:4 with lines title "OpenFOAM-deal.II", \
	     "watchpoints/openfoam-fenics.log" using 1:4 with lines title "OpenFOAM-FEniCS", \
	     "watchpoints/su2-calculix.log" using 1:4 with lines title "SU2-CalculiX", \
	     "watchpoints/su2-dealii.log" using 1:4 with lines title "SU2-deal.II", \
	     "watchpoints/su2-fenics.log" using 1:4 with lines title "SU2-FEniCS", \
	     "watchpoints/nutils-calculix.log" using 1:4 with lines title "Nutils-CalculiX", \
	     "watchpoints/nutils-dealii.log" using 1:4 with lines title "Nutils-deal.II", \
	     "watchpoints/nutils-fenics.log" using 1:4 with lines title "Nutils-FEniCS"
EOF
