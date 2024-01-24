#!/bin/sh
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'displacement at the middle of the tube                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'displacement [m]'
	set term pngcairo enhanced size 900,654
	set output "images/tutorials-elastic-tube-3d-displacement-all-watchpoints.png"
	plot "$1" using 1:5 with lines title "OpenFOAM-CalculiX circumferential", \
	     "$1" using 1:6 with lines title "OpenFOAM-CalculiX radial", \
	     "$1" using 1:7 with lines title "OpenFOAM-CalculiX axial"
EOF
