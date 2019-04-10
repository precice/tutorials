There are 2 main things missing in the structural FEniCS-FEniCS solver missing:

1. Correct expressions describing Dirichlet and Neumann boundary
    * For Dirichlet BC we use displacement
    * Should we use Force or Stress for Neumann BC?
2. Loop for solving the problem. It is a static problem, which is generally not handled by preCICE, but can be circumvented by either reformulation of the problem into:
		
		laplacian(u) - f = du/dt,

or by approaching the problem in an iterative way.
