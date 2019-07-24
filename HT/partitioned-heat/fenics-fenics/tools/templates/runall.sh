python3 ../../../../heat.py -d -wr {{ wr_dirichlet }} {{ wr_neumann }} -dT {{ window_size }} -cpl {{ coupling_scheme }} -g {{gamma}} -tol {{error_tolerance}} &
python3 ../../../../heat.py -n -wr {{ wr_dirichlet }} {{ wr_neumann }} -dT {{ window_size }} -cpl {{ coupling_scheme }} -g {{gamma}} -tol {{error_tolerance}} 
