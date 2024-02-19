#!/bin/sh

gnuplot -p << EOF
  set grid
  set term pngcairo enhanced size 900,654
  set output "images/tutorials-elastic-tube-1d-all.png"

  set multiplot layout 2, 1
  set title 'Diameter of elastic-tube at x=5'

  set ylabel 'diameter[m]'
  plot \
    "fluid-cpp/precice-Fluid-watchpoint-Middle.log" using 1:4 with linespoints title "cpp", \
    "fluid-python/precice-Fluid-watchpoint-Middle.log" using 1:4 with linespoints lw 2 title "python", \
    "fluid-rust/precice-Fluid-watchpoint-Middle.log" using 1:4 with linespoints title "rust"

  set title 'Pressure of elastic-tube at x=5'
  set ylabel 'pressure'
  plot \
    "fluid-cpp/precice-Fluid-watchpoint-Middle.log" using 1:5 with linespoints title "cpp", \
    "fluid-python/precice-Fluid-watchpoint-Middle.log" using 1:5 with linespoints title "python", \
    "fluid-rust/precice-Fluid-watchpoint-Middle.log" using 1:5 with linespoints title "rust"

	set xlabel 'time [s]'

EOF
