#! /bin/bash
gnuplot -p << EOF
set grid
set title 'Displacement of the Flap Tip'
set xlabel 'Time [s]'
set ylabel 'X-Displacement [m]'
plot "point1.watchpoint.txt" using 1:8 with lines title ""
EOF
