#!/bin/sh

if [ "${1:-}" = "" ]; then
    echo "No target directory specified. Please specify the directory of the solid participant containing the watchpoint, e.g. ./plot-displacements.sh solid-calculix."
    exit 1
fi

FILE="$1/precice-Solid-watchpoint-Tube-Midpoint.log"

if [ ! -f "$FILE" ]; then
    echo "Unable to locate the watchpoint file (precice-Solid-watchpoint-Tube-Midpoint.log) in the specified solid directory '${1}'. Make sure the specified directory matches the solid participant you used for the calculations."
    exit 1
fi
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'displacement of the tube midpoint'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'displacement of the tube midpoint [m]'
	set term pngcairo enhanced size 900,654
	set output "images/tutorials-elastic-tube-3d-displacement-watchpoints.png"
	plot "$1/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:5 with lines title "$1 tube midpoint circumferential", \
	     "$1/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:6 with lines title "$1 tube midpoint radial", \
	     "$1/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:7 with lines title "$1 tube midpoint axial"
EOF

echo "Plot saved in images/tutorials-elastic-tube-3d-displacement-watchpoints.png"
