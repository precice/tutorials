#!/bin/sh

gnuplot -p << EOF
  set grid
  set term pngcairo enhanced size 900,654
  set output "images/tutorials-elastic-tube-1d-all.png"

  set multiplot layout 2, 1
  set title 'Diameter of elastic-tube at x=5'

  set ylabel 'diameter[m]'
  plot "runs/solid-cpp-fluid-rust/middle.log" using 1:4 with lines title "cpp-rust", \
       "runs/solid-cpp-fluid-python/middle.log" using 1:4 with lines title "cpp-python", \
       "runs/solid-cpp-fluid-cpp/middle.log" using 1:4 with lines title "cpp-cpp", \
       "runs/solid-python-fluid-rust/middle.log" using 1:4 with lines title "python-rust", \
       "runs/solid-python-fluid-python/middle.log" using 1:4 with lines title "python-python", \
       "runs/solid-python-fluid-cpp/middle.log" using 1:4 with lines title "python-cpp", \
       "runs/solid-rust-fluid-rust/middle.log" using 1:4 with lines title "rust-rust", \
       "runs/solid-rust-fluid-python/middle.log" using 1:4 with lines title "rust-python", \
       "runs/solid-rust-fluid-cpp/middle.log" using 1:4 with lines title "rust-cpp"

  set title 'Pressure of elastic-tube at x=5'
  set ylabel 'pressure'
  plot "runs/solid-cpp-fluid-rust/middle.log" using 1:5 with lines title "cpp-rust", \
       "runs/solid-cpp-fluid-python/middle.log" using 1:5 with lines title "cpp-python", \
       "runs/solid-cpp-fluid-cpp/middle.log" using 1:5 with lines title "cpp-cpp", \
       "runs/solid-python-fluid-rust/middle.log" using 1:5 with lines title "python-rust", \
       "runs/solid-python-fluid-python/middle.log" using 1:5 with lines title "python-python", \
       "runs/solid-python-fluid-cpp/middle.log" using 1:5 with lines title "python-cpp", \
       "runs/solid-rust-fluid-rust/middle.log" using 1:5 with lines title "rust-rust", \
       "runs/solid-rust-fluid-python/middle.log" using 1:5 with lines title "rust-python", \
       "runs/solid-rust-fluid-cpp/middle.log" using 1:5 with lines title "rust-cpp"

	set xlabel 'time [s]'

EOF
