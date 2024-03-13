#!/bin/sh
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'x-displacement of the flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'x-displacement [m]'
	set term pngcairo enhanced size 900,654
	set output "images/tutorials-perpendicular-flap-displacement-all-watchpoints.png"
	plot "watchpoints/fake-fenics.log" using 1:4 with lines title "fake-FEniCS", \
	     "watchpoints/fake-calculix.log" using 1:4 with lines title "fake-CalculiX"
EOF
