#!/usr/bin/env sh                                                                    

if [ "${1:-}" = "" ]; then
    echo "No target directory specified. Please specify the directory of the solid participant containing the watchpoint, e.g. ./plot-displacement.sh solid-dealii."
    exit 1
fi

FILE="$1/precice-Solid-watchpoint-Flap-Tip.log"

if [ ! -f "$FILE" ]; then
    echo "Unable to locate the watchpoint file (precice-Solid-watchpoint-Flap-Tip.log) in the specified solid directory '${1}'. Make sure the specified directory matches the solid participant you used for the calculations."
    exit 1
fi

gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'y-displacement of the flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'y-displacement [m]'                                                 
	set linestyle  1 lt 2 lc 1 # red-dashed                                         
	plot "$1/precice-Solid-watchpoint-Flap-Tip.log" using 1:5 with lines title "$1"
EOF
