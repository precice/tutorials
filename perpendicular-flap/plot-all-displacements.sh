#!/bin/sh
 
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
