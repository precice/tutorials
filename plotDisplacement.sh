#! /bin/bash
gnuplot -p << EOF
set grid
set title 'Displacement of the Tip of the Flap'
set xlabel 'Time'
set ylabel 'Displacement'
plot "point1.watchpoint.txt" using 1:5 with line lc 'red' title ""
EOF
