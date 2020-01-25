#! /bin/bash
gnuplot -p << EOF
set grid
set title 'Displacement of the Flap Tip'
set xlabel 'Time [s]'
set ylabel 'Y-Displacement [m]'
set xrange [34:35]
set yrange [-0.15:0.2]
plot "referenceY-disp" title "Turek-Hron FSI benchmark FSI2" lc rgb "black", "precice-Calculix-watchpoint-point1.log" using 1:9 with lines title "present work" lc rgb "red" lw 2
EOF

