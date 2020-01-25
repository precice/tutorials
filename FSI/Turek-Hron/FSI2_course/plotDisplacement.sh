#! /bin/bash
gnuplot -p << EOF
set grid
set title 'Displacement of the Flap Tip'
set xlabel 'Time [s]'
set ylabel 'Y-Displacement [m]'
plot "precice-Calculix-watchpoint-point1.log" using 1:9 with lines title "present work" lc rgb "red" lw 2
EOF

