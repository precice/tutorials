import numpy as np
from mpi4py import MPI
import basix.ufl
from petsc4py import PETSc
import ufl
from dolfinx import fem, io, mesh as msh, default_scalar_type
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, create_vector, set_bc, LinearProblem
from dolfinx.mesh import create_rectangle
import basix
from fenicsxprecice import Adapter, CouplingMesh


# geometry
nx = 100
ny = 25
nz = 1

y_top = 0
y_bottom = y_top - .25
x_left = 0
x_right = x_left + 1

fenics_dt = 0.01  # time step size


def top_boundary(x):
    tol = 1E-14
    return np.isclose(x[1], y_top, tol)


def bottom_boundary(x):
    tol = 1E-14
    return np.isclose(x[1], y_bottom, tol)


class initial_value():
    def __init__(self, constant):
        self.constant = constant

    def __call__(self, x):
        return np.full(x[0].shape, self.constant)


class GradientSolver:
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    The solver has been changed since the original version from the link above introduces larger errors

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

    def compute(self, u, k):
        L = fem.form(ufl.inner(-k * ufl.grad(u), self.v) * ufl.dx)
        b = create_vector(fem.extract_function_spaces(L))
        assemble_vector(b, L)
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        self.solver.solve(b, self.returnValue.x.petsc_vec)
        return self.returnValue


p0 = (x_left, y_bottom)
p1 = (x_right, y_top)

mesh = create_rectangle(MPI.COMM_WORLD, [np.asarray(p0), np.asarray(p1)], [nx, ny], msh.CellType.triangle)
V = fem.functionspace(mesh, ('P', 2))
# for the vector function space
element = basix.ufl.element("Lagrange", mesh.topology.cell_name(), 1, shape=(mesh.geometry.dim,))
V_g = fem.functionspace(mesh, element)
W, map_to_W = V_g.sub(1).collapse()

gradient_solver = GradientSolver(mesh, V_g)

alpha = 1  # m^2/s, https://en.wikipedia.org/wiki/Thermal_diffusivity
k = 100  # kg * m / s^3 / K, https://en.wikipedia.org/wiki/Thermal_conductivity

# We will only exchange flux in y direction on coupling interface. No initialization necessary.
flux_y = fem.Function(W)

# Define initial value
u_n = fem.Function(V)
u_n.name = "T"
u_n.interpolate(initial_value(310))

tdim = mesh.topology.dim
fdim = tdim - 1
mesh.topology.create_connectivity(fdim, tdim)
dofs_coupling = fem.locate_dofs_geometrical(V, top_boundary)
dofs_bottom = fem.locate_dofs_geometrical(V, bottom_boundary)

# Adapter definition and initialization
precice = Adapter(adapter_config_filename="precice-adapter-config.json", mpi_comm=MPI.COMM_WORLD)

# top_boundary is coupling boundary
coupling_mesh = CouplingMesh("Solid-Mesh", top_boundary, {"Temperature": V}, {"Heat-Flux": flux_y})
precice.initialize([coupling_mesh])

# boundary function for the coupling interface
coupling_function = fem.Function(V)

# Assigning appropriate dt
precice_dt = precice.get_max_time_step_size()
dt = np.min([fenics_dt, precice_dt])

# Define variational problem
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
F = u * v / dt * ufl.dx + alpha * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx - u_n * v / dt * ufl.dx

# apply constant Dirichlet boundary condition at bottom edge
# apply Dirichlet boundary condition on coupling interface
bcs = [fem.dirichletbc(coupling_function, dofs_coupling), fem.dirichletbc(default_scalar_type(310), dofs_bottom, V)]

a = fem.form(ufl.lhs(F))
L = fem.form(ufl.rhs(F))

A = assemble_matrix(a, bcs=bcs)
A.assemble()
b = create_vector(fem.extract_function_spaces(L))
uh = fem.Function(V)
solver = PETSc.KSP().create(mesh.comm)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)


# Time-stepping
t = 0

vtxwriter = io.VTXWriter(MPI.COMM_WORLD, f"output_solid.bp", [u_n])
vtxwriter.write(t)

n = 0

flux = fem.Function(V_g)

while precice.is_coupling_ongoing():

    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, 0)

    precice_dt = precice.get_max_time_step_size()
    dt = np.min([fenics_dt, precice_dt])
    precice.read_data(coupling_mesh.get_name(), "Temperature", dt, coupling_function)

    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
    assemble_vector(b, L)

    apply_lifting(b, [a], [bcs])
    set_bc(b, bcs)

    # Solve linear problem
    solver.solve(b, uh.x.petsc_vec)

    # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
    flux = gradient_solver.compute(u_n, k)
    flux_y.interpolate(flux.sub(1))
    precice.write_data(coupling_mesh.get_name(), "Heat-Flux", flux_y)

    precice.advance(dt)

    if precice.requires_reading_checkpoint():
        u_cp, t_cp, _ = precice.retrieve_checkpoint()
        u_n.x.array[:] = u_cp.x.array
        t = t_cp
    else:  # update solution
        # Update solution at previous time step (u_n)
        u_n.x.array[:] = uh.x.array
        t += float(dt)
        n += 1
        if n % 20 == 0:
            vtxwriter.write(t)


precice.finalize()
vtxwriter.close()
