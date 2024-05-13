#!/usr/bin/env sh

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

WATCHPOINTS_DIR="./reference-results/watchpoints/"

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (selected combinations)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-selected-watchpoints.png"
    plot "${WATCHPOINTS_DIR}/openfoam-calculix-v2404.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
         "${WATCHPOINTS_DIR}/openfoam-dealii-v2404.log" using 1:4 with lines title "OpenFOAM-deal.II", \
         "${WATCHPOINTS_DIR}/openfoam-fenics-v2404.log" using 1:4 with lines title "OpenFOAM-FEniCS", \
         "${WATCHPOINTS_DIR}/su2-calculix-v2404.log" using 1:4 with lines title "SU2-CalculiX", \
         "${WATCHPOINTS_DIR}/su2-dealii-v2404.log" using 1:4 with lines title "SU2-deal.II", \
         "${WATCHPOINTS_DIR}/su2-fenics-v2404.log" using 1:4 with lines title "SU2-FEniCS", \
         "${WATCHPOINTS_DIR}/nutils-calculix-v2404.log" using 1:4 with lines title "Nutils-CalculiX", \
         "${WATCHPOINTS_DIR}/nutils-dealii-v2404.log" using 1:4 with lines title "Nutils-deal.II"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (comparison of flow solvers)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-flow-comparison-watchpoints.png"
    plot "${WATCHPOINTS_DIR}/openfoam-calculix-v2404.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
         "${WATCHPOINTS_DIR}/su2-calculix-v2404.log" using 1:4 with lines title "SU2-CalculiX", \
         "${WATCHPOINTS_DIR}/nutils-calculix-v2404.log" using 1:4 with lines title "Nutils-CalculiX", \
         "${WATCHPOINTS_DIR}/fake-calculix-v2404.log" using 1:4 with lines title "fake-CalculiX"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (incompressible flow with OpenFOAM)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-openfoam-watchpoints.png"
    plot "${WATCHPOINTS_DIR}/openfoam-calculix-v2404.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
         "${WATCHPOINTS_DIR}/openfoam-dealii-v2404.log" using 1:4 with lines title "OpenFOAM-deal.II", \
         "${WATCHPOINTS_DIR}/openfoam-dune-v2404.log" using 1:4 with lines title "OpenFOAM-DUNE", \
         "${WATCHPOINTS_DIR}/openfoam-fenics-v2404.log" using 1:4 with lines title "OpenFOAM-FEniCS", \
         "${WATCHPOINTS_DIR}/openfoam-nutils-v2404.log" using 1:4 with lines title "OpenFOAM-Nutils", \
         "${WATCHPOINTS_DIR}/openfoam-solids4foam-v2404.log" using 1:4 with lines title "OpenFOAM-solids4Foam"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (compressible flow with SU2)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-su2-watchpoints.png"
    plot "${WATCHPOINTS_DIR}/su2-calculix-v2404.log" using 1:4 with lines title "SU2-CalculiX", \
         "${WATCHPOINTS_DIR}/su2-dealii-v2404.log" using 1:4 with lines title "SU2-deal.II", \
         "${WATCHPOINTS_DIR}/su2-dune-v2404.log" using 1:4 with lines title "SU2-DUNE", \
         "${WATCHPOINTS_DIR}/su2-fenics-v2404.log" using 1:4 with lines title "SU2-FEniCS", \
         "${WATCHPOINTS_DIR}/su2-nutils-v2404.log" using 1:4 with lines title "SU2-nutils", \
         "${WATCHPOINTS_DIR}/su2-solids4foam-v2404.log" using 1:4 with lines title "SU2-solids4Foam"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (incompressible flow with Nutils)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-nutils-watchpoints.png"
    plot "${WATCHPOINTS_DIR}/nutils-calculix-v2404.log" using 1:4 with lines title "Nutils-CalculiX", \
         "${WATCHPOINTS_DIR}/nutils-dealii-v2404.log" using 1:4 with lines title "Nutils-deal.II", \
         "${WATCHPOINTS_DIR}/nutils-dune-v2404.log" using 1:4 with lines title "Nutils-DUNE", \
         "${WATCHPOINTS_DIR}/nutils-fenics-v2404.log" using 1:4 with lines title "Nutils-FEniCS", \
         "${WATCHPOINTS_DIR}/nutils-nutils-v2404.log" using 1:4 with lines title "Nutils-Nutils", \
         "${WATCHPOINTS_DIR}/nutils-solids4foam-v2404.log" using 1:4 with lines title "Nutils-solids4Foam"
EOF

gnuplot -p << EOF
    set grid
    set title 'x-displacement of the flap tip (dummy force data with fluid-fake)'
    set xlabel 'time [s]'
    set ylabel 'x-displacement [m]'
    set term pngcairo enhanced size 900,654
    set output "images/tutorials-perpendicular-flap-displacement-fake-watchpoints.png"
    plot "${WATCHPOINTS_DIR}/fake-calculix-v2404.log" using 1:4 with lines title "fake-CalculiX", \
         "${WATCHPOINTS_DIR}/fake-dealii-v2404.log" using 1:4 with lines title "fake-deal.II", \
         "${WATCHPOINTS_DIR}/fake-dune-v2404.log" using 1:4 with lines title "fake-DUNE", \
         "${WATCHPOINTS_DIR}/fake-fenics-v2404.log" using 1:4 with lines title "fake-FEniCS", \
         "${WATCHPOINTS_DIR}/fake-nutils-v2404.log" using 1:4 with lines title "fake-Nutils", \
         "${WATCHPOINTS_DIR}/fake-solids4foam-v2404.log" using 1:4 with lines title "fake-solids4Foam"
EOF
