#! /bin/bash                                                                    
gnuplot -p << EOF                                                               
        set grid                                                                        
        set title 'Displacement of the rigid body tip'                                        
        set xlabel 'time [s]'                                                           
        set ylabel 'y-displacement [m]'                                                 
        plot "solid-cpp/precice-Solid-watchpoint-Flap-Tip.log" using 1:5 with lines notitle
EOF
