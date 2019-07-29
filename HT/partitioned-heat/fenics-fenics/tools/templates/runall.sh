python3 ../../../../heat.py -d {{domain_decomposition_dirichlet}} -wr {{ wr_dirichlet }} {{ wr_neumann }} -dT {{ window_size }} -cpl {{ coupling_scheme }} -g {{gamma}} -tol {{error_tolerance}} &
python3 ../../../../heat.py -n {{domain_decomposition_neumann}} -wr {{ wr_dirichlet }} {{ wr_neumann }} -dT {{ window_size }} -cpl {{ coupling_scheme }} -g {{gamma}} -tol {{error_tolerance}}
