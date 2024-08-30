#!/bin/sh

if [ "${1:-}" = "" ]; then
    echo "No target directory specified. Please specify the directory of the participant containing the watchpoint, e.g. ./plot-displacement.sh coil-python."
    exit 1
fi

FILE="$1/precice-Capacitor-watchpoint-VoltageCurrent.log"

if [ ! -f "$FILE" ]; then
    echo "Unable to locate the watchpoint file (precice-Capacitor-watchpoint-VoltageCurrent.log) in the specified participant directory '${1}'. Make sure the specified directory matches the participant you used for the calculations."
    exit 1
fi

gnuplot -p << EOF
	set grid
	set title 'Voltage and current'
	set xlabel 'time [s]'
	set ylabel 'Voltage / Current'
	plot "$1/precice-Capacitor-watchpoint-VoltageCurrent.log" using 1:4 with linespoints title "Voltage", "$1/precice-Capacitor-watchpoint-VoltageCurrent.log" using 1:5 with linespoints title "Current"
EOF
