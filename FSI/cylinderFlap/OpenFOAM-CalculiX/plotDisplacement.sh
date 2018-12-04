#! /bin/bash
gnuplot -p << EOF
set grid
set title 'Displacement of the Flap Tip'
set xlabel 'Time [s]'
set ylabel 'Y-Displacement [m]'
plot "point1.watchpoint.txt" using 1:9 with lines title ""
EOF

