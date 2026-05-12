"""
The basic example is taken from "Langtangen, Hans Petter, and Anders Logg. Solving PDEs in Python: The FEniCS
Tutorial I. Springer International Publishing, 2016."

The example code has been extended with preCICE API calls and mixed boundary conditions to allow for a Dirichlet-Neumann
coupling of two separate heat equations. It also has been adapted to be compatible with FEniCSx.

The original source code can be found on https://jsdokken.com/dolfinx-tutorial/chapter2/heat_equation.html

Heat equation with Dirichlet conditions. (Dirichlet problem)
  u'= Laplace(u) + f  in the unit square [0,1] x [0,1]
  u = u_C             on the coupling boundary at x = 1
  u = u_D             on the remaining boundary
  u = u_0             at t = 0
  u = 1 + x^2 + alpha*y^2 + \beta*t
  f = beta - 2 - 2*alpha

Heat equation with mixed boundary conditions. (Neumann problem)
  u'= Laplace(u) + f  in the shifted unit square [1,2] x [0,1]
  du/dn = f_N         on the coupling boundary at x = 1
  u = u_D             on the remaining boundary
  u = u_0             at t = 0
  u = 1 + x^2 + alpha*y^2 + \beta*t
  f = beta - 2 - 2*alpha
"""
import argparse
import numpy as np
from mpi4py import MPI
import basix.ufl
from petsc4py import PETSc
import ufl
from dolfinx import fem, io, mesh as msh, geometry
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, create_vector, set_bc
import basix
from fenicsxprecice import Adapter, CouplingMesh
from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart
from problem_setup import get_geometry


def interface(x):
    return np.isclose(x[0], 1)


class GradientSolver:
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    The given legacy FEniCS code has been modified to be compatible with FEniCSx.
    Additionally, the GradientSolver reuses the system matrix and the solver

    :param V_g: Vector function space
    :param u: solution where gradient is to be determined
    """

    def __init__(self, domain, V_g):
        self.domain = domain,
        self.V_g = V_g

        w = ufl.TrialFunction(V_g)
        self.v = ufl.TestFunction(V_g)
        a = fem.form(ufl.inner(w, self.v) * ufl.dx)
        self.A = assemble_matrix(a)
        self.A.assemble()

        self.solver = PETSc.KSP().create(domain.comm)
        self.solver.setOperators(self.A)
        self.solver.setType(PETSc.KSP.Type.PREONLY)
        self.solver.getPC().setType(PETSc.PC.Type.LU)

        self.returnValue = fem.Function(V_g)

    def compute(self, u):
        L = fem.form(ufl.inner(ufl.grad(u), self.v) * ufl.dx)
        b = create_vector(fem.extract_function_spaces(L))
        assemble_vector(b, L)
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        self.solver.solve(b, self.returnValue.x.petsc_vec)
        return self.returnValue


# Parse arguments
parser = argparse.ArgumentParser(description="Solving heat equation for simple or complex interface case")
parser.add_argument("participantName", help="Name of the solver.", type=str, choices=[p.value for p in ProblemType])
parser.add_argument("-e", "--error-tol", help="set error tolerance", type=float, default=10**-9,)
args = parser.parse_args()
participant_name = args.participantName
error_tol = args.error_tol

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

t = 0
fenics_dt = 0.1
alpha = 3
beta = 1.2

# define the domain
if participant_name == ProblemType.DIRICHLET.value:
    problem = ProblemType.DIRICHLET
    domain_part = DomainPart.LEFT
elif participant_name == ProblemType.NEUMANN.value:
    problem = ProblemType.NEUMANN
    domain_part = DomainPart.RIGHT

# create domain and function space
domain, coupling_boundary, remaining_boundary = get_geometry(domain_part, comm)
V = fem.functionspace(domain, ("Lagrange", 2))
element = basix.ufl.element("Lagrange", domain.topology.cell_name(), 1, shape=(domain.geometry.dim,))
V_g = fem.functionspace(domain, element)
W, map_to_W = V_g.sub(0).collapse()
V_coup = None
if problem is ProblemType.DIRICHLET:
    V_coup = V
else:
    V_coup = W

# Define the exact solution


class exact_solution():
    def __init__(self, alpha, beta, t):
        self.alpha = alpha
        self.beta = beta
        self.t = t

    def __call__(self, x):
        return 1 + x[0]**2 + self.alpha * x[1]**2 + self.beta * self.t


u_exact = exact_solution(alpha, beta, t)

gradient_solver = GradientSolver(domain, V_g)

# Define the boundary condition
bcs = []
u_D = fem.Function(V)
u_D.interpolate(u_exact)
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
# dofs for the coupling boundary
dofs_boundary = fem.locate_dofs_geometrical(V, coupling_boundary)
# dofs for the remaining boundary. Can be directly set to u_D
dofs_remaining = fem.locate_dofs_geometrical(V, remaining_boundary)
bc_D = fem.dirichletbc(u_D, dofs_remaining)
bcs.append(bc_D)

# dofs (and the corresponding coordinates) for the coupling boundary (coupling function space)
# this is later used to read data from preCICE and update boundary condition
dofs_coupling = fem.locate_dofs_geometrical(V_coup, coupling_boundary)
dofs_coupling_coordinates = V_coup.tabulate_dof_coordinates()[dofs_coupling]

if problem is ProblemType.DIRICHLET:
    # Define flux in x direction
    f_N = fem.Function(W)
    f_N.interpolate(lambda x: 2 * x[0])

u_n = fem.Function(V)  # IV and solution u for the n-th time step
u_n.interpolate(u_exact)

precice, precice_dt, initial_data = None, 0.0, None
if problem is ProblemType.DIRICHLET:
    precice = Adapter(adapter_config_filename="precice-adapter-config-D.json", mpi_comm=comm)
else:
    precice = Adapter(adapter_config_filename="precice-adapter-config-N.json", mpi_comm=comm)

coupling_mesh = None
if problem is ProblemType.DIRICHLET:
    coupling_mesh = CouplingMesh("Dirichlet-Mesh", coupling_boundary, {"Temperature": V}, {"Heat-Flux": f_N})
    precice.initialize([coupling_mesh])
elif problem is ProblemType.NEUMANN:
    coupling_mesh = CouplingMesh("Neumann-Mesh", coupling_boundary, {"Heat-Flux": V}, {"Temperature": u_D})
    precice.initialize([coupling_mesh])

# get precice's dt
precice_dt = precice.get_max_time_step_size()
dt = np.min([fenics_dt, precice_dt])

# Define the variational formualation
# As $f$ is a constant independent of $t$, we can define it as a constant.
f = fem.Constant(domain, beta - 2 - 2 * alpha)
# We can now create our variational formulation, with the bilinear form `a` and  linear form `L`.
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
F = u * v * ufl.dx + dt * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx - (u_n + dt * f) * v * ufl.dx


name_read = None
if problem is ProblemType.DIRICHLET:
    name_read = "Temperature"
    read_mesh = "Dirichlet-Mesh"
else:
    name_read = "Heat-Flux"
    read_mesh = "Neumann-Mesh"

coupling_function = fem.Function(V_coup)

if problem is ProblemType.DIRICHLET:
    # modify Dirichlet boundary condition on coupling interface
    bc_coup = fem.dirichletbc(coupling_function, dofs_boundary)
    bcs.append(bc_coup)
if problem is ProblemType.NEUMANN:
    # modify Neumann boundary condition on coupling interface, modify weak
    # form correspondingly
    F += dt * coupling_function * v * ufl.ds

a = fem.form(ufl.lhs(F))
L = fem.form(ufl.rhs(F))

# ## Create the matrix and vector for the linear problem
# To ensure that we are solving the variational problem efficiently, we
# will create several structures which can reuse data, such as matrix
# sparisty patterns. Especially note as the bilinear form `a` is
# independent of time, we only need to assemble the matrix once.
A = assemble_matrix(a, bcs=bcs)
A.assemble()
b = create_vector(fem.extract_function_spaces(L))
uh = fem.Function(V)

# ## Define a linear variational solver
# We will use [PETSc](https://www.mcs.anl.gov/petsc/) to solve the
# resulting linear algebra problem. We use the Python-API `petsc4py` to
# define the solver. We will use a linear solver.
solver = PETSc.KSP().create(domain.comm)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)

if problem is ProblemType.DIRICHLET:
    flux = fem.Function(V_g)

# values of the boundary conditions are set to the end of the timestep
u_exact.t += dt
u_D.interpolate(u_exact)

bb_tree = geometry.bb_tree(domain, domain.geometry.dim)
cell_candidates = geometry.compute_collisions_points(bb_tree, dofs_coupling_coordinates)
cell_candidates = geometry.compute_colliding_cells(domain, cell_candidates, dofs_coupling_coordinates)
cc = np.array(list(set(cell_candidates.array)))

f_err = fem.Function(V)
vtxwriter = io.VTXWriter(MPI.COMM_WORLD, f"output_{problem.name}.bp", [f_err])
vtxwriter.write(t)

while precice.is_coupling_ongoing():

    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, 0)

    precice_dt = precice.get_max_time_step_size()
    dt = np.min([fenics_dt, precice_dt])

    # Update the coupling expression with the new read data
    precice.read_data(read_mesh, name_read, dt, coupling_function)

    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
    assemble_vector(b, L)

    # Apply Dirichlet boundary condition to the vector (according to the tutorial, the lifting operation is used to preserve the symmetry of the matrix)
    # Boundary condition bc should be updated by u_D.interpolate above, since
    # this function is wrapped into the bc object
    apply_lifting(b, [a], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)

    set_bc(b, bcs)

    # Solve linear problem
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()

    # Write data to preCICE according to which problem is being solved
    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem reads temperature and writes flux on boundary to Neumann problem
        flux = gradient_solver.compute(uh)
        flux_x = fem.Function(W)
        flux_x.interpolate(flux.sub(0))
        flux_x.x.scatter_forward()
        precice.write_data(coupling_mesh.get_name(), "Heat-Flux", flux_x)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem reads flux and writes temperature on boundary to Dirichlet problem
        precice.write_data(coupling_mesh.get_name(), "Temperature", uh)

    precice.advance(dt)
    precice_dt = precice.get_max_time_step_size()

    # revert to checkpoint
    if precice.requires_reading_checkpoint():
        u_cp, t_cp, _ = precice.retrieve_checkpoint()
        u_n.x.array[:] = u_cp.x.array
        t = t_cp
    else:  # update solution
        # Update solution at previous time step (u_n)
        u_n.x.array[:] = uh.x.array
        f_err.x.array[:] = np.abs(u_n.x.array - u_D.x.array)
        t += float(dt)
        vtxwriter.write(t)

    if precice.is_time_window_complete():
        u_ref = fem.Function(V)
        u_ref.interpolate(u_D)
        error, error_pointwise = compute_errors(u_n, u_ref, total_error_tol=error_tol)
        print("t = %.2f: L2 error on domain = %.3g" % (t, error))
        # Update Dirichlet BC
        u_exact.t += dt
        u_D.interpolate(u_exact)
        # TODO: update time dependent f (as soon as it is time dependent)!

precice.finalize()

vtxwriter.close()
