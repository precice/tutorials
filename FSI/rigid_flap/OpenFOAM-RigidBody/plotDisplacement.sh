#! /bin/bash                                                                    
gnuplot -p << EOF                                                               
set grid                                                                        
set title 'Displacement of the rigid body tip'                                        
set xlabel 'time [s]'                                                           
set ylabel 'y-displacement [m]'                                                 
plot "precice-Solid-watchpoint-flap-tip.log" using 1:5 notitle with lines 
EOF

