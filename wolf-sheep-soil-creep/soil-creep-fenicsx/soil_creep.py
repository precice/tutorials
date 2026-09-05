"""
Linear diffusion equation of soil with Dirichlet boundary.
  u' = div(D*grad(u))    in the square [0,19] x [0,19]
  u = 0                  on the bottom boundary at y=0
"""

import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import ufl
import numpy as np
from petsc4py import PETSc
from mpi4py import MPI
from dolfinx import default_scalar_type, fem, mesh
from dolfinx.fem.petsc import (
    assemble_vector,
    assemble_matrix,
    create_vector,
    apply_lifting,
    set_bc,
)
from fenicsxprecice import Adapter, CouplingMesh

def bottom_boundary(x):
    return np.isclose(x[1], 0)

# Volume coupling
def coupling_boundary(x):
    return np.ones(x.shape[1], dtype=bool)


initial_soil_depth = 0.3

# Set the soil-thickness scale for limiting creep where little soil is available
hstar = 0.2

# Set parameters for two soil creep coefficients: slow (full grass cover) and fast (partial or "eaten" grass cover)
fast_creep = 0.1
slow_creep = 0.001

# Set grid size (same as W-S-G model's grid)
# Adjust refinement by changing number of nodes ny/nx (e.g. double refinement by setting ny = nx = 40).
nx = ny = 20
length_x = length_y = 19.0
dx = length_x / (nx - 1)

fenics_dt = 0.2 * dx * dx / fast_creep

comm = MPI.COMM_WORLD

# Create domain and function space
domain = mesh.create_rectangle(
    comm,
    [np.array([0.0, 0.0]), np.array([length_x, length_y])],
    [nx - 1, ny - 1],
    mesh.CellType.triangle,
)
V = fem.functionspace(domain, ("Lagrange", 1))
Q = fem.functionspace(domain, ("DG", 0))

# u is elev function
u_n = fem.Function(V)
u_n.interpolate(lambda x: 0.1 * x[1])

initial_elev = fem.Function(V)
initial_elev.x.array[:] = u_n.x.array

# D is diffusivity function / creep coefficient
D = fem.Function(Q)

# Functions for coupling
soil = fem.Function(V)
soil.x.array[:] = initial_soil_depth

grass = fem.Function(V)

# Auxiliary functions
grass_diffusivity = fem.Function(V)
soil_scaling = fem.Function(V)

# Dirichlet boundary condition
bottom_dofs = fem.locate_dofs_geometrical(V, bottom_boundary)
bc = fem.dirichletbc(default_scalar_type(0), bottom_dofs, V)

# preCICE setup
participant = Adapter(adapter_config_filename="precice-adapter-config.json", mpi_comm=comm)

coupling_mesh = CouplingMesh("Soil-Creep-Mesh", coupling_boundary, {"Grass": V}, {"Soil": soil})
participant.initialize([coupling_mesh])

# Define the variational formualation
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
dt_constant = fem.Constant(domain, default_scalar_type(fenics_dt))

a = u * v * ufl.dx + dt_constant * D * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = u_n * v * ufl.dx

bilinear_form = fem.form(a)
linear_form = fem.form(L)

# Create the matrix and vector for the linear problem
b = create_vector(fem.extract_function_spaces(linear_form))
uh = fem.Function(V)

# Define a linear variational solver
solver = PETSc.KSP().create(domain.comm)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)


while participant.is_coupling_ongoing():
    precice_dt = participant.get_max_time_step_size()
    dt = np.min([fenics_dt, precice_dt])
    dt_constant.value = default_scalar_type(dt)

    participant.read_data("Soil-Creep-Mesh", "Grass", dt, grass)

    grass_diffusivity.x.array[:] = np.where(
        grass.x.array == 1,
        fast_creep,
        slow_creep,
    )
    grass_diffusivity.x.scatter_forward()

    soil_scaling.x.array[:] = 1.0 - np.exp(-soil.x.array / hstar)
    soil_scaling.x.scatter_forward()

    D_expr = fem.Expression(
        grass_diffusivity * soil_scaling,
        Q.element.interpolation_points,
    )
    D.interpolate(D_expr)

    # Reassemble problem because D is updated
    A = assemble_matrix(bilinear_form, bcs=[bc])
    A.assemble()
    solver.setOperators(A)

    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
    assemble_vector(b, linear_form)

    apply_lifting(b, [bilinear_form], [[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)

    set_bc(b, [bc])

    # Solve linear problem
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()

    soil.x.array[:] += uh.x.array - u_n.x.array
    soil.x.scatter_forward()

    # Update solution at previous time step (u_n)
    u_n.x.array[:] = uh.x.array

    participant.write_data("Soil-Creep-Mesh", "Soil", soil)

    participant.advance(dt)


participant.finalize()



dofs_coupling_coordinates = V.tabulate_dof_coordinates()[:, :2].copy()

def values_to_grid(values):
    grid = np.full((ny, nx), np.nan)

    for coord, value in zip(dofs_coupling_coordinates, values):
        i = int(round(coord[0] / length_x * (nx - 1)))
        j = int(round(coord[1] / length_y * (ny - 1)))
        grid[j, i] = value

    return grid


# Mask for closed boundaries (mimic Landlab)
closed = np.zeros((ny, nx), dtype=bool)

closed[1:-1, 0] = True
closed[1:-1, -1] = True
closed[-1, :] = True

# Calculate and plot the erosion/deposition patterns
ero_dep = u_n.x.array - initial_elev.x.array
ero_dep_grid = values_to_grid(ero_dep)

maxchange = np.max(np.abs(ero_dep_grid))

cmap = mpl.colormaps["coolwarm_r"].copy()
cmap.set_bad("black")

plt.figure()
plt.imshow(
    np.ma.array(ero_dep_grid, mask=closed),
    origin="lower",
    vmin=-maxchange,
    vmax=maxchange,
    cmap=cmap
)
plt.colorbar(label="Depth of soil accumulation (+) or loss (-), m")
plt.savefig("output/erosion_deposition_patterns.png")
plt.close()

# Soil thickness
soil_grid = values_to_grid(soil.x.array)

cmap = mpl.colormaps["pink"].copy()
cmap.set_bad("black")

plt.figure()
plt.imshow(
    np.ma.array(soil_grid, mask=closed),
    origin="lower",
    cmap=cmap,
)
plt.colorbar(label="Soil thickness, m")
plt.savefig("output/soil_thickness.png")
plt.close()

# Ground cover
gm_grid = values_to_grid(grass.x.array)

cmap = mpl.colormaps["YlGn"].copy()
cmap.set_bad("black")

plt.figure()
plt.imshow(
    np.ma.array(gm_grid, mask=closed),
    origin="lower",
    cmap=cmap,
    vmin=1,
    vmax=2,
)
plt.colorbar(label="Ground cover (1 = bare, 2 = grass)")
plt.savefig("output/grass_map.png")
plt.close()
