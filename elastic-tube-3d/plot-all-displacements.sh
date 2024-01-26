#!/bin/sh
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'displacement at the middle of the tube'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'displacement [m]'
	set term pngcairo enhanced size 900,654
	set output "images/tutorials-elastic-tube-3d-displacement-all-watchpoints.png"
	plot "solid-calculix/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:5 with lines title "OpenFOAM-CalculiX circumferential", \
	     "solid-calculix/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:6 with lines title "OpenFOAM-CalculiX radial", \
	     "solid-calculix/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:7 with lines title "OpenFOAM-CalculiX axial", \
	     "solid-fenics/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:5 with lines title "OpenFOAM-FEniCS circumferential", \
	     "solid-fenics/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:6 with lines title "OpenFOAM-FEniCS radial", \
	     "solid-fenics/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:7 with lines title "OpenFOAM-FEniCS axial
EOF
