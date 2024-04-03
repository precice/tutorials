#!/bin/sh

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (solid-openfoam)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-openfoam-convergence-watchpoints.png"
    plot "watchpoints/openfoam-openfoam-v3.1.0.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 6x15", \
         "watchpoints/openfoam-openfoam-6-100.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 6x100", \
         "watchpoints/openfoam-openfoam-6-150.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 6x150", \
         "watchpoints/openfoam-openfoam-10-100.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 10x100", \
         "watchpoints/openfoam-openfoam-10-150.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 10x150", \
         "watchpoints/openfoam-solids4foam-v3.1.0.log" using 1:4 with lines title "OpenFOAM-solids4Foam"
EOF

        #  "watchpoints/openfoam-openfoam-6-60.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 6x60", \
        #  "watchpoints/openfoam-openfoam-12-30.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 12x30", \
        #  "watchpoints/openfoam-openfoam-24-60.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 24x60", \
        #  "watchpoints/openfoam-openfoam-20-100.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 20x100", \
        #  "watchpoints/openfoam-openfoam-20-150.log" using 1:4 with lines title "OpenFOAM-OpenFOAM 20x150", \
