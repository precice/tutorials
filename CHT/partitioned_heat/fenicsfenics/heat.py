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


class Tests(Enum):
    SIN = 1
    FENICS = 2


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


def fluxes_from_temperature_full_domain(F, bcs, mesh, dt, hy):
    V = FunctionSpace(mesh, 'CG', 1)
    fluxes_vector = assemble(F)
    fluxes = Function(V)
    fluxes.vector()[:] = - fluxes_vector[:] * dt / hy
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

# Create mesh and define function space
nx = ny = 20

T = 1.0  # final time
num_steps = 10  # number of time steps
dt = T / num_steps  # time step size
alpha = 3  # parameter alpha
beta = 1.3  # parameter beta
x_coupling = .7  # x coordinate of coupling interface
y_bottom = 0
y_top = 1
hy = (y_top - y_bottom) / (ny)

TEST_NAME = Tests.SIN

if TEST_NAME is Tests.FENICS:
    lam = 1
elif TEST_NAME is Tests.SIN:
    lam = .01

lam_c = Constant(lam)

if problem is ProblemType.DIRICHLET:
    p0 = Point(0, y_bottom)
    p1 = Point(x_coupling, y_top)
elif problem is ProblemType.NEUMANN:
    p0 = Point(x_coupling, y_bottom)
    p1 = Point(2, y_top)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
if TEST_NAME is Tests.FENICS:
    u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=dt)
elif TEST_NAME is Tests.SIN:
    u_D = Expression('exp(-5 * t *lam * pi * pi / 4) * (500 * sin(pi/2 * x[0]) * sin (pi * x[1]))', t=0, lam=lam,
                     degree=2)
u_D_function = interpolate(u_D, V)
# Define flux on coupling interface (grad(u_D) in normal direction)
if TEST_NAME is Tests.FENICS:
    f_N = Expression('2 * x[0]', degree=1)
elif TEST_NAME is Tests.SIN:
    f_N = Expression('500 * pi/2 * cos(pi/2 * x[0]) * sin (pi * x[1])', degree=1)
f_N_function = interpolate(f_N, V)

coupling_boundary = CouplingBoundary()
remaining_boundary = ComplementaryBoundary(coupling_boundary)

# todo put all the function calls below into one?
coupling = fenics_adapter.Coupling(config_file_name, solver_name)
coupling.set_coupling_mesh(mesh, coupling_boundary, coupling_mesh_name)
if problem is ProblemType.DIRICHLET:
    coupling.set_read_field(read_data_name, u_D_function)
    coupling.set_write_field(write_data_name, f_N_function)
elif problem is ProblemType.NEUMANN:
    coupling.set_read_field(read_data_name, f_N_function)
    coupling.set_write_field(write_data_name, u_D_function)
coupling.initialize_data()
# todo until here.

bcs = [DirichletBC(V, u_D, remaining_boundary)]
# Define initial value
u_n = interpolate(u_D, V)
# u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
if TEST_NAME is Tests.FENICS:
    f = Constant(beta - 2 - 2 * alpha)
elif TEST_NAME is Tests.SIN:
    f = Constant(0)

F = u * v * dx + dt * lam_c * dot(grad(u), grad(v)) * dx - (u_n + dt * f) * v * dx

if problem is ProblemType.DIRICHLET:
    # apply Dirichlet boundary condition on coupling interface
    bcs.append(coupling.create_coupling_dirichlet_boundary_condition(V))
if problem is ProblemType.NEUMANN:
    # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
    F += coupling.create_coupling_neumann_boundary_condition(v)

a, L = lhs(F), rhs(F)

"""
normal = dolfin.FacetNormal(mesh)
mesh_function = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
mesh_function.set_all(0)
remaining_boundary.mark(mesh_function, 1)
non_coupling_ds = dolfin.Measure('ds', domain=mesh, subdomain_data=mesh_function)
"""

# Time-stepping
u = Function(V)
F_alternative = (u - (u_n + dt * f)) / dt * v * dx + lam_c * dot(grad(u), grad(v)) * dx #- dot(normal, grad(u)) * v * non_coupling_ds
u.rename("Temperature", "")
t = coupling.precice_tau
u_D.t = t
assert (dt == coupling.precice_tau)

file_out = File("out/%s.pvd" % solver_name)
ref_out = File("out/ref%s.pvd" % solver_name)

while coupling.is_coupling_ongoing():

    # Compute solution
    solve(a == L, u, bcs)

    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
        fluxes = fluxes_from_temperature_full_domain(F_alternative, bcs, mesh, dt, hy)
        coupling.exchange_data(fluxes, dt)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem obtains sends temperature on boundary to Dirichlet problem
        coupling.exchange_data(u, dt)

    is_converged = coupling.check_convergence()

    if is_converged:
        # plot solution
        # plot(u)
        # plt.pause(.1)
        # Compute error at vertices
        u_e = interpolate(u_D, V)
        u_e.rename("reference", " ")
        # error = np.abs(u_e.vector().get_local() - u.vector().get_local()).max()
        error = assemble(inner(u_e - u, u_e - u) * dx)
        print('t = %.2f: error = %.3g' % (t, error))
        # Update previous solution
        file_out << u
        ref_out << u_e
        # Update current time
        t += coupling.precice_tau
        # Update dirichlet BC
        u_D.t = t
        u_n.assign(u)

# Hold plot
coupling.finalize()
