#!/bin/sh

if [ "${1:-}" = "" ]; then
    echo "No target directory specified. Please specify the directory of the solid participant containing the watchpoint, e.g. ./plot-displacement.sh solid-calculix."
    exit 1
fi

FILE="$1/precice-Solid-watchpoint-flaptip.log"

if [ ! -f "$FILE" ]; then
    echo "Unable to locate the watchpoint file (precice-Solid-watchpoint-flaptip.log) in the specified solid directory '${1}'. Make sure the specified directory matches the solid participant you used for the calculations."
    exit 1
fi
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'x-displacement of the flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'x-displacement [m]'                                                 
	plot "$1/precice-Solid-watchpoint-flaptip.log" using 1:9 with lines title "$1"
EOF
