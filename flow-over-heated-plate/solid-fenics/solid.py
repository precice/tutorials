"""
Solid plate participant in flow-over-plate tutorial using FEniCS
"""

from __future__ import print_function, division
from fenics import Function, SubDomain, RectangleMesh, BoxMesh, FunctionSpace, VectorFunctionSpace, Point, \
    Expression, Constant, DirichletBC, \
    TrialFunction, TestFunction, File, solve, plot, lhs, rhs, grad, inner, dot, dx, ds, interpolate, project, \
    near, MeshFunction, MPI
from fenicsprecice import Adapter
import numpy as np


class ComplementaryBoundary(SubDomain):
    """Determines if a point is at the complementary boundary with tolerance of
    1E-14.
    :func inside(): returns True if point belongs to the boundary, otherwise
                    returns False
    """

    def __init__(self, subdomain):
        self.complement = subdomain
        SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and not self.complement.inside(x, on_boundary):
            return True
        else:
            return False


class TopBoundary(SubDomain):
    """Determines if the point is at the top boundary with tolerance of 1E-14.
    :func inside(): returns True if point belongs to the boundary, otherwise
                    returns False
    """

    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and near(x[1], y_top, tol):
            return True
        else:
            return False


class BottomBoundary(SubDomain):
    """Determines if the point is at the bottom boundary with tolerance of
    1E-14.

    :func inside(): returns True if point belongs to the boundary, otherwise
                    returns False
    """

    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and near(x[1], y_bottom, tol):
            return True
        else:
            return False


def determine_heat_flux(V_g, u, k, flux):
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    :param V_g: Vector function space
    :param u: solution where gradient is to be determined
    :param k: thermal conductivity
    :param flux: returns calculated flux into this value
    """

    w = TrialFunction(V_g)
    v = TestFunction(V_g)

    a = inner(w, v) * dx
    L = inner(-k * grad(u), v) * dx
    solve(a == L, flux)


# Create mesh and define function space
nx = 100
ny = 25
nz = 1

fenics_dt = 0.01  # time step size
dt_out = 0.2  # interval for writing VTK files
y_top = 0
y_bottom = y_top - .25
x_left = 0
x_right = x_left + 1

p0 = Point(x_left, y_bottom, 0)
p1 = Point(x_right, y_top, 1)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 1)
V_g = VectorFunctionSpace(mesh, 'P', 1)

alpha = 1  # m^2/s, https://en.wikipedia.org/wiki/Thermal_diffusivity
k = 100  # kg * m / s^3 / K, https://en.wikipedia.org/wiki/Thermal_conductivity

# Define boundary condition
u_D = Constant('310')
u_D_function = interpolate(u_D, V)
# We will only exchange flux in y direction on coupling interface. No initialization necessary.
V_flux_y = V_g.sub(1)

coupling_boundary = TopBoundary()
bottom_boundary = BottomBoundary()

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("T", "")

# Adapter definition and initialization
precice = Adapter(adapter_config_filename="precice-adapter-config.json")

precice_dt = precice.initialize(coupling_boundary, read_function_space=V, write_object=V_flux_y)

# Create a FEniCS Expression to define and control the coupling boundary values
coupling_expression = precice.create_coupling_expression()

# Assigning appropriate dt
dt = Constant(0)
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
F = u * v / dt * dx + alpha * dot(grad(u), grad(v)) * dx - u_n * v / dt * dx

# apply constant Dirichlet boundary condition at bottom edge
# apply Dirichlet boundary condition on coupling interface
bcs = [DirichletBC(V, coupling_expression, coupling_boundary), DirichletBC(V, u_D, bottom_boundary)]

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
t = 0
u_D.t = t + dt

# mark mesh w.r.t ranks
ranks = File("Solid/VTK/ranks%s.pvd.pvd" % precice.get_participant_name())
mesh_rank = MeshFunction("size_t", mesh, mesh.topology().dim())
mesh_rank.set_all(MPI.rank(MPI.comm_world))
mesh_rank.rename("myRank", "")
ranks << mesh_rank

# Create output file
file_out = File("Solid/VTK/%s.pvd" % precice.get_participant_name())
file_out << u_n

print("output vtk for time = {}".format(float(t)))
n = 0

fluxes = Function(V_g)
fluxes.rename("Fluxes", "")

while precice.is_coupling_ongoing():

    if precice.is_action_required(precice.action_write_iteration_checkpoint()):  # write checkpoint
        precice.store_checkpoint(u_n, t, n)

    read_data = precice.read_data()

    # Update the coupling expression with the new read data
    precice.update_coupling_expression(coupling_expression, read_data)

    dt.assign(np.min([fenics_dt, precice_dt]))

    # Compute solution
    solve(a == L, u_np1, bcs)

    # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
    determine_heat_flux(V_g, u_np1, k, fluxes)
    fluxes_y = fluxes.sub(1)  # only exchange y component of flux.
    precice.write_data(fluxes_y)

    precice_dt = precice.advance(dt(0))

    if precice.is_action_required(precice.action_read_iteration_checkpoint()):  # roll back to checkpoint
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:  # update solution
        u_n.assign(u_np1)
        t += float(dt)
        n += 1

    if precice.is_time_window_complete():
        tol = 10e-5  # we need some tolerance, since otherwise output might be skipped.
        if abs((t + tol) % dt_out) < 2 * tol:  # output if t is a multiple of dt_out
            print("output vtk for time = {}".format(float(t)))
            file_out << u_n

    # Update dirichlet BC
    u_D.t = t + float(dt)

# Hold plot
precice.finalize()
