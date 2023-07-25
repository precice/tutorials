import sys
import itertools
import numpy as np

rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

# Time stepping
dt_3D = 1e-1            # time step of 3D mechanics
end_time = 30.0         # end time of the simulation 
output_interval = dt_3D # time interval between outputs

# Material parameters
pmax = 7.3                                                  # maximum active stress
rho = 10                                                    # density of the muscle
material_parameters = [3.176e-10, 1.813, 1.075e-2, 1.0]     # [c1, c2, b, d]

# Meshes
ex_x, ex_y, ex_z = 3.0, 3.0, 12.0               # extent of muscle
el_x, el_y, el_z = 3, 3, 12                     # number of elements
bs_x, bs_y, bs_z = 2*el_x+1, 2*el_y+1, 2*el_z+1 # quadratic basis functions

fb_x, fb_y = 10, 10         # number of fibers
fb_points = 100             # number of points per fiber
fiber_direction = [0, 0, 1] # direction of fiber in element

meshes = { # create 3D mechanics mesh
    "mesh3D": {
        "nElements":            [el_x, el_y, el_z],
        "physicalExtent":       [ex_x, ex_y, ex_z],
        "physicalOffset":       [0, 0, 0],
        "logKey":               "mesh3D",
        "inputMeshIsGlobal":    True,
        "nRanks":               n_ranks
    }
}

# Boundary conditions
dirichlet_bc = {} # fix z=0 with dirichlet boundary conditions
for x in range(bs_x):
    for y in range(bs_y):
        dirichlet_bc[x + y*bs_x] = [0.0, 0.0, 0.0, None, None, None]

neumann_bc = [] # add pulling force to z=el_z with neumann boundary conditions
neumann_force = 0
for x in range(el_x):
    for y in range(el_y):
        neumann_bc += [{
            "element": x + y*el_x + (el_z-1)*el_y*el_x, 
            "constantVector": [0, 0, neumann_force], 
            "face": "2+"
        }]