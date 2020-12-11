#! /bin/bash                                                                    
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'y-displacement of the flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'y-displacement [m]'                                                 
	set linestyle  1 lt 2 lc 1 # red-dashed                                         
	plot "precice-Solid-watchpoint-flap_tip.log" using 1:5 with lines
EOF
