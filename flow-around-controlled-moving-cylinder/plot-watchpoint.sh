#!/bin/sh

if [ "${1:-}" = "" ]; then
    echo "No target directory specified. Please specify the directory of the Solid participant containing the watchpoint, e.g. ./plot-watchpoint.sh solid-python."
    exit 1
fi

FILE="$1/precice-Solid-watchpoint-Spring.log"

if [ ! -f "$FILE" ]; then
	echo "Unable to locate the watchpoint file (*.log) in the specified directory '${1}'. Make sure the specified directory matches the fluid or solid participant you used for the calculations."
	exit 1
fi
 
gnuplot -p <<EOF                                                         
	set grid                                                                        
	set title 'Displacement spring'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'delta y [m]'                                                 
	plot "$FILE" using 1:5 with lines title "$1"
EOF

gnuplot -p <<EOF                                             
	set grid                                                                        
	set title 'Displacement cylinder'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'delta y [m]'                                                 
	plot "$FILE" using 1:7 with lines title "$1"
EOF

gnuplot -p <<EOF                                                             
	set grid                                                                        
	set title 'Drag Force'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'force y [N]'                                                 
	plot "$FILE" using 1:8 with lines title "$1"
EOF

gnuplot -p <<EOF                                                             
	set grid                                                                        
	set title 'Lift Force'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'force y [N]'                                                 
	plot "$FILE" using 1:9 with lines title "$1"
EOF
