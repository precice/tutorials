#!/bin/sh                                                                    
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'y-displacement of the flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'y-displacement [m]'                                                 
	set linestyle  1 lt 2 lc 1 # red-dashed                                         
	plot "solid-dealii/precice-Solid-watchpoint-Flap-Tip.log" using 1:5 with lines notitle
EOF
