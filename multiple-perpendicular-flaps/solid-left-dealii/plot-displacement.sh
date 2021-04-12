#!/bin/sh
gnuplot -p << EOF
set grid
set title 'Displacement of the Flap Tip'
set xlabel 'Time [s]'
set ylabel 'X-Displacement [m]'
plot "precice-Solid1-watchpoint-flap_tip.log" using 1:4 title 'Top displacemement' with lines
EOF

