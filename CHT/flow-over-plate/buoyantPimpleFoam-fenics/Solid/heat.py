"""
The basic example is taken from "Langtangen, Hans Petter, and Anders Logg. Solving PDEs in Python: The FEniCS
Tutorial I. Springer International Publishing, 2016."

The example code has been extended with preCICE API calls and mixed boundary conditions to allow for a Dirichlet-Neumann
coupling of two separate heat equations.

The original source code can be found on https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft03_heat.py.

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

from __future__ import print_function, division
from fenics import Function, SubDomain, RectangleMesh, BoxMesh, FunctionSpace, Point, Expression, Constant, DirichletBC, \
    TrialFunction, TestFunction, File, solve, plot, lhs, rhs, grad, inner, dot, dx, ds, assemble, interpolate, project, near
from fenicsadapter import Adapter
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



def fluxes_from_temperature_full_domain(F, V, k):
    """Computes flux from weak form (see p.3 in Toselli, Andrea, and Olof
    Widlund. Domain decomposition methods-algorithms and theory. Vol. 34.
    Springer Science & Business Media, 2006.).

    :param F: weak form with known u^{n+1}
    :param V: function space
    :param k: thermal conductivity
    :return: fluxes function
    """
    fluxes_vector = assemble(F)  # assemble weak form -> evaluate integral
    v = TestFunction(V)
    fluxes = Function(V)  # create function for flux
    area = assemble(v * ds).get_local()
    for i in range(area.shape[0]):
        if area[i] != 0:  # put weight from assemble on function
            fluxes.vector()[i] = - k * fluxes_vector[i] / area[i]  # scale by surface area
        else:
            assert(abs(fluxes_vector[i]) < 10**-10)  # for non surface parts, we expect zero flux
            fluxes.vector()[i] = - k * fluxes_vector[i]
    return fluxes


# Create mesh and define function space
nx = 100
ny = 25
nz = 1

fenics_dt = 0.01  # time step size
dt_out = 0.2
y_top = 0
y_bottom = y_top - .25
x_left = 0
x_right = x_left + 1

p0 = Point(x_left, y_bottom, 0)
p1 = Point(x_right, y_top, 1)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 1)

alpha = 1  # m^2/s, https://en.wikipedia.org/wiki/Thermal_diffusivity
k = 100  # kg * m / s^3 / K, https://en.wikipedia.org/wiki/Thermal_conductivity

# Define boundary condition
u_D = Constant('310')
u_D_function = interpolate(u_D, V)
# Define flux in x direction on coupling interface (grad(u_D) in normal direction)
f_N = Constant('0')
f_N_function = interpolate(f_N, V)

coupling_boundary = TopBoundary()
bottom_boundary = BottomBoundary()

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("T", "")

#start preCICE adapter
precice = Adapter()
precice_dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=u_D_function, write_field=f_N_function,
                   u_n=u_n)
dt = np.min([fenics_dt, precice_dt])  # todo we could also consider deciding on time stepping size inside the adapter
bcs = [DirichletBC(V, u_D, bottom_boundary)]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
F = u * v / dt * dx + alpha * dot(grad(u), grad(v)) * dx - u_n * v / dt * dx

# apply Dirichlet boundary condition on coupling interface
bcs.append(precice.create_coupling_dirichlet_boundary_condition(V))

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
F_known_u = u_np1 * v / dt * dx + alpha * dot(grad(u_np1), grad(v)) * dx - u_n * v / dt * dx
t = 0
u_D.t = t + dt

file_out = File("Solid/VTK/%s.pvd" % precice._solver_name)
n=0

while precice.is_coupling_ongoing():

    # Compute solution
    solve(a == L, u_np1, bcs)

    # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
    fluxes = fluxes_from_temperature_full_domain(F_known_u, V, k)
    t, n, precice_timestep_is_complete, precice_dt = precice.advance(fluxes, u_np1, u_n, t, dt, n)

    dt = np.min([fenics_dt, precice_dt])  # todo we could also consider deciding on time stepping size inside the adapter

    if precice_timestep_is_complete:
        if abs(t % dt_out) < 10e-5 or abs(t % dt_out) < 10e-5:  # just a very complicated way to only produce output if t is a multiple of dt_out
            file_out << u_n

    # Update dirichlet BC
    u_D.t = t + dt

file_out << u_n

# Hold plot
precice.finalize()
