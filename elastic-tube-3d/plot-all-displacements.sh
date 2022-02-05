#!/bin/sh
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'radial displacement of the tube midpoint'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'radial displacement [m]'
	set term pngcairo enhanced size 900,654
	set output "images/tutorials-elastic-tube-3d-displacement-all-watchpoints.png"
	plot "solid-calculix/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
	     "solid-fenics/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:4 with lines title "OpenFOAM-FEniCS"
EOF
