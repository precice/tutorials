#! /bin/bash
gnuplot -p << EOF
set grid
set title 'Displacement of the Flap Tip'
set xlabel 'Time [s]'
set ylabel 'X-Displacement [m]'
plot "precice-Calculix-watchpoint-point1.log" using 1:5 smooth cumulative title ""
EOF
