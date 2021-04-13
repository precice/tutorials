#!/bin/sh
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'x-displacement of the flap tip'                                        
	set xlabel 'time [s]'                                                           
	set ylabel 'x-displacement [m]'                                                 
	plot "openfoam-calculix.log" using 1:4 with lines title "OpenFOAM-CalculiX", \
	     "openfoam-dealii.log" using 1:4 with lines title "OpenFOAM-deal.II", \
	     "openfoam-fenics.log" using 1:4 with lines title "OpenFOAM-FEniCS", \
	     "su2-calculix.log" using 1:4 with lines title "SU2-CalculiX", \
	     "su2-dealii.log" using 1:4 with lines title "SU2-deal.II", \
	     "su2-fenics.log" using 1:4 with lines title "SU2-FEniCS", \
	     "nutils-dealii.log" using 1:4 with lines title "Nutils-deal.II"
EOF
