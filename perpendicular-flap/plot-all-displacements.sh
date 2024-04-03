#!/bin/sh

# This script cannot be used as-is and is meant to generate the picture
# images/tutorials-perpendicular-flap-displacement-all-watchpoints.png
# rendered on https://precice.org/tutorials-perpendicular-flap.html.
#
# It plots watchpoints of specific combinations of Fluid and Solid participants,
# stored in `images/`. To use this script:
# 1. For each combination you want to plot:
#    1. Run the tutorial with that combination
#    2. Copy the respective watchpoint file to watchpoints/, with names <fluid>-<solid>.log
#    3. Clean the tutorial
# 2. Edit the script to plot the files you want, with the corresponding titles,
#    adding one line per combination.
# 3. Call the script with ./plot-all-displacements.sh, from this directory.
#
# If you are only interested in a subset of combinations, remove the respective lines.

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-all-watchpoints.png"
    plot "watchpoints/openfoam-calculix-v3.1.0.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
         "watchpoints/openfoam-dealii-v3.1.0.log" using 1:4 with lines title "OpenFOAM-deal.II", \
         "watchpoints/openfoam-dune-v3.1.0.log" using 1:4 with lines title "OpenFOAM-DUNE", \
         "watchpoints/openfoam-fenics-v3.1.0.log" using 1:4 with lines title "OpenFOAM-FEniCS", \
         "watchpoints/openfoam-nutils-v3.1.0.log" using 1:4 with lines title "OpenFOAM-Nutils", \
         "watchpoints/openfoam-openfoam-v3.1.0.log" using 1:4 with lines title "OpenFOAM-OpenFOAM", \
         "watchpoints/openfoam-solids4foam-v3.1.0.log" using 1:4 with lines title "OpenFOAM-solids4Foam", \
         "watchpoints/su2-calculix-v3.1.0.log" using 1:4 with lines title "SU2-CalculiX", \
         "watchpoints/su2-dealii-v3.1.0.log" using 1:4 with lines title "SU2-deal.II", \
         "watchpoints/su2-dune-v3.1.0.log" using 1:4 with lines title "SU2-DUNE", \
         "watchpoints/su2-fenics-v3.1.0.log" using 1:4 with lines title "SU2-FEniCS", \
         "watchpoints/su2-nutils-v3.1.0.log" using 1:4 with lines title "SU2-nutils", \
         "watchpoints/su2-openfoam-v3.1.0.log" using 1:4 with lines title "SU2-OpenFOAM", \
         "watchpoints/su2-solids4foam-v3.1.0.log" using 1:4 with lines title "SU2-solids4Foam", \
         "watchpoints/nutils-calculix-v3.0.0.log" using 1:4 with lines title "Nutils-CalculiX", \
         "watchpoints/nutils-dealii-v3.1.0.log" using 1:4 with lines title "Nutils-deal.II", \
         "watchpoints/nutils-nutils-v3.0.0.log" using 1:4 with lines title "Nutils-Nutils"
EOF


gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (selected combinations)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-selected-watchpoints.png"
    plot "watchpoints/openfoam-calculix-v3.1.0.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
         "watchpoints/openfoam-dealii-v3.1.0.log" using 1:4 with lines title "OpenFOAM-deal.II", \
         "watchpoints/openfoam-fenics-v3.1.0.log" using 1:4 with lines title "OpenFOAM-FEniCS", \
         "watchpoints/su2-calculix-v3.1.0.log" using 1:4 with lines title "SU2-CalculiX", \
         "watchpoints/su2-dealii-v3.1.0.log" using 1:4 with lines title "SU2-deal.II", \
         "watchpoints/su2-fenics-v3.1.0.log" using 1:4 with lines title "SU2-FEniCS", \
         "watchpoints/nutils-calculix-v3.0.0.log" using 1:4 with lines title "Nutils-CalculiX", \
         "watchpoints/nutils-dealii-v3.1.0.log" using 1:4 with lines title "Nutils-deal.II"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (comparison of flow solvers)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-flow-comparison-watchpoints.png"
    plot "watchpoints/openfoam-calculix-v3.1.0.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
         "watchpoints/su2-calculix-v3.1.0.log" using 1:4 with lines title "SU2-CalculiX", \
         "watchpoints/nutils-calculix-v3.0.0.log" using 1:4 with lines title "Nutils-CalculiX", \
         "watchpoints/fake-calculix-v3.1.0.log" using 1:4 with lines title "fake-CalculiX"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (incompressible flow with OpenFOAM)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-openfoam-watchpoints.png"
    plot "watchpoints/openfoam-calculix-v3.1.0.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
         "watchpoints/openfoam-dealii-v3.1.0.log" using 1:4 with lines title "OpenFOAM-deal.II", \
         "watchpoints/openfoam-dune-v3.1.0.log" using 1:4 with lines title "OpenFOAM-DUNE", \
         "watchpoints/openfoam-fenics-v3.1.0.log" using 1:4 with lines title "OpenFOAM-FEniCS", \
         "watchpoints/openfoam-nutils-v3.1.0.log" using 1:4 with lines title "OpenFOAM-Nutils", \
         "watchpoints/openfoam-solids4foam-v3.1.0.log" using 1:4 with lines title "OpenFOAM-solids4Foam"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (compressible flow with SU2)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-su2-watchpoints.png"
    plot "watchpoints/su2-calculix-v3.1.0.log" using 1:4 with lines title "SU2-CalculiX", \
         "watchpoints/su2-dealii-v3.1.0.log" using 1:4 with lines title "SU2-deal.II", \
         "watchpoints/su2-dune-v3.1.0.log" using 1:4 with lines title "SU2-DUNE", \
         "watchpoints/su2-fenics-v3.1.0.log" using 1:4 with lines title "SU2-FEniCS", \
         "watchpoints/su2-nutils-v3.1.0.log" using 1:4 with lines title "SU2-nutils", \
         "watchpoints/su2-solids4foam-v3.1.0.log" using 1:4 with lines title "SU2-solids4Foam"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (incompressible flow with Nutils)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-nutils-watchpoints.png"
    plot "watchpoints/nutils-calculix-v3.0.0.log" using 1:4 with lines title "Nutils-CalculiX", \
         "watchpoints/nutils-dealii-v3.1.0.log" using 1:4 with lines title "Nutils-deal.II", \
         "watchpoints/nutils-nutils-v3.0.0.log" using 1:4 with lines title "Nutils-Nutils"
EOF

# Not currently included in the Nutils plots due to long simulation time:
    #      "watchpoints/nutils-fenics.log" using 1:4 with lines title "Nutils-FEniCS", \
    #      "watchpoints/nutils-dune.log" using 1:4 with lines title "Nutils-DUNE", \
    #      "watchpoints/nutils-openfoam.log" using 1:4 with lines title "Nutils-OpenFOAM", \
    #      "watchpoints/nutils-solids4foam.log" using 1:4 with lines title "Nutils-solids4Foam", \

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (dummy force data with fluid-fake)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-fake-watchpoints.png"
    plot "watchpoints/fake-calculix-v3.1.0.log" using 1:4 with lines title "fake-CalculiX", \
         "watchpoints/fake-dealii-v3.1.0.log" using 1:4 with lines title "fake-deal.II", \
         "watchpoints/fake-dune-v3.1.0.log" using 1:4 with lines title "fake-DUNE", \
         "watchpoints/fake-fenics-v3.1.0.log" using 1:4 with lines title "fake-FEniCS", \
         "watchpoints/fake-nutils-v3.1.0.log" using 1:4 with lines title "fake-Nutils", \
         "watchpoints/fake-solids4foam-v3.1.0.log" using 1:4 with lines title "fake-solids4Foam"
EOF
