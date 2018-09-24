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
from fenics import Function, SubDomain, RectangleMesh, FunctionSpace, Point, Expression, Constant, DirichletBC, \
    TrialFunction, TestFunction, File, solve, plot, lhs, rhs, grad, inner, dot, dx, ds, assemble, interpolate, project, near
import dolfin
import numpy as np
from enum import Enum
from matplotlib import pyplot as plt
import fenics_adapter
import argparse


class ProblemType(Enum):
    """
    Enum defines problem type. Details see above.
    """
    DIRICHLET = 1  # Dirichlet problem
    NEUMANN = 2  # Neumann problem


class ComplementaryBoundary(SubDomain):
    def __init__(self, subdomain):
        self.complement = subdomain
        SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and not self.complement.inside(x, on_boundary):
            return True
        else:
            return False


class CouplingBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and near(x[0], x_coupling, tol) and not near(x[1], y_bottom, tol) and not near(x[1], y_top, tol):
            return True
        else:
            return False


def fluxes_from_temperature_full_domain(F, V, hy):
    """
    compute flux from weak form (see p.3 in Toselli, Andrea, and Olof Widlund. Domain decomposition methods-algorithms and theory. Vol. 34. Springer Science & Business Media, 2006.)
    :param F: weak form with known u^{n+1}
    :param V: function space
    :param hy: spatial resolution perpendicular to flux direction
    :return:
    """
    fluxes_vector = assemble(F)  # assemble weak form -> evaluate integral
    fluxes = Function(V)  # create function for flux
    fluxes.vector()[:] = fluxes_vector[:] / hy  # put weight from assemble on function, scale function by spatial resolution
    return fluxes


parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml precice configuration file.", type=str)
parser.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest='dirichlet', action='store_true')
parser.add_argument("-n", "--neumann", help="create a neumann problem", dest='neumann', action='store_true')

try:
    args = parser.parse_args()
except SystemExit:
    raise Exception("No config file name given. Did you forget adding the precice configuration file as an argument?")

config_file_name = args.configurationFileName

# coupling parameters
if args.dirichlet:
    problem = ProblemType.DIRICHLET
if args.neumann:
    problem = ProblemType.NEUMANN
if args.dirichlet and args.neumann:
    raise Exception("you can only choose either a dirichlet problem (option -d) or a neumann problem (option -n)")
if not (args.dirichlet or args.neumann):
    raise Exception("you have to choose either a dirichlet problem (option -d) or a neumann problem (option -n)")

# Create mesh and define function space
nx = ny = 30

if problem is ProblemType.DIRICHLET:
    solver_name = "HeatDirichlet"
    coupling_mesh_name = "DirichletNodes"
    read_data_name = "Temperature"
    write_data_name = "Flux"
elif problem is ProblemType.NEUMANN:
    solver_name = "HeatNeumann"
    coupling_mesh_name = "NeumannNodes"
    read_data_name = "Flux"
    write_data_name = "Temperature"

T = 1.0  # final time
num_steps = 10  # number of time steps
dt = T / num_steps  # time step size
alpha = 3  # parameter alpha
beta = 1.3  # parameter beta
y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = .7  # x coordinate of coupling interface
hy = (y_top - y_bottom) / ny

if problem is ProblemType.DIRICHLET:
    p0 = Point(x_left, y_bottom)
    p1 = Point(x_coupling, y_top)
elif problem is ProblemType.NEUMANN:
    p0 = Point(x_coupling, y_bottom)
    p1 = Point(x_right, y_top)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=dt)
u_D_function = interpolate(u_D, V)
# Define flux in x direction on coupling interface (grad(u_D) in normal direction)
f_N = Expression('2 * x[0]', degree=1)
f_N_function = interpolate(f_N, V)

coupling_boundary = CouplingBoundary()
remaining_boundary = ComplementaryBoundary(coupling_boundary)

coupling = fenics_adapter.Coupling(config_file_name, solver_name)
coupling.set_coupling_mesh(mesh, coupling_boundary, coupling_mesh_name)
if problem is ProblemType.DIRICHLET:
    coupling.set_read_field(read_data_name, u_D_function)
    coupling.set_write_field(write_data_name, f_N_function)
elif problem is ProblemType.NEUMANN:
    coupling.set_read_field(read_data_name, f_N_function)
    coupling.set_write_field(write_data_name, u_D_function)
coupling.initialize_data()

bcs = [DirichletBC(V, u_D, remaining_boundary)]
# Define initial value
u_n = interpolate(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2 * alpha)
F = u * v * dx + dt * dot(grad(u), grad(v)) * dx - (u_n + dt * f) * v * dx

if problem is ProblemType.DIRICHLET:
    # apply Dirichlet boundary condition on coupling interface
    bcs.append(coupling.create_coupling_dirichlet_boundary_condition(V))
if problem is ProblemType.NEUMANN:
    # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
    F += coupling.create_coupling_neumann_boundary_condition(v)

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
F_known_u = u_np1 * v * dx + dt * dot(grad(u_np1), grad(v)) * dx - (u_n + dt * f) * v * dx
u_np1.rename("Temperature", "")
t = 0
u_D.t = t + coupling.precice_tau
assert (dt == coupling.precice_tau)

file_out = File("out/%s.pvd" % solver_name)
ref_out = File("out/ref%s.pvd" % solver_name)

while coupling.is_coupling_ongoing():

    # Compute solution
    solve(a == L, u_np1, bcs)

    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
        fluxes = fluxes_from_temperature_full_domain(F_known_u, V, hy)
        coupling.exchange_data(fluxes, dt)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem obtains sends temperature on boundary to Dirichlet problem
        coupling.exchange_data(u_np1, dt)

    is_converged = coupling.check_convergence()

    if is_converged:
        # Compute error at vertices
        u_e = interpolate(u_D, V)
        u_e.rename("reference", " ")
        error = assemble(inner(u_e - u_np1, u_e - u_np1)/(u_e * u_e) * dx)
        assert (error < 10e-4)
        print('t = %.2f: error = %.3g' % (t, error))
        # Update previous solution
        file_out << u_np1
        ref_out << u_e
        # Update current time
        t += coupling.precice_tau
        # Update dirichlet BC
        u_D.t = t + coupling.precice_tau
        u_n.assign(u_np1)

# Hold plot
coupling.finalize()
