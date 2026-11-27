#!/usr/bin/env sh

if [ "${1:-}" = "" ]; then
    echo "No target directory specified. Please specify the directory of the solid participant containing the watchpoint, e.g. ./plot-displacement.sh solid-openradioss."
    exit 1
fi

FILE_CENTER="$1/precice-Solid-watchpoint-Center-Front-Tip.log"

if [ ! -f "$FILE_CENTER" ]; then
    echo "Unable to locate the watchpoint file ${FILE_CENTER}. Make sure the specified directory matches the solid participant you used for the calculations."
    exit 1
fi
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'x-displacement of the central flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'x-displacement [m]'                                                 
	plot "$1/precice-Solid-watchpoint-Center-Front-Tip.log" using 1:5 with lines title "$1"
EOF

FILE_LEFT="$1/precice-Solid-watchpoint-Left-Back-Tip.log"

if [ ! -f "$FILE_LEFT" ]; then
    echo "Unable to locate the watchpoint file ${FILE_LEFT}. Make sure the specified directory matches the solid participant you used for the calculations."
    exit 1
fi
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'x-displacement of the left flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'x-displacement [m]'                                                 
	plot "$1/precice-Solid-watchpoint-Left-Back-Tip.log" using 1:5 with lines title "$1"
EOF