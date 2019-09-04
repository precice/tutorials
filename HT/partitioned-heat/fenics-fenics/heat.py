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
    TrialFunction, TestFunction, File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate, assemble, project, near, VectorFunctionSpace
from enum import Enum
from fenicsadapter import Adapter, ExactInterpolationExpression, GeneralInterpolationExpression
from errorcomputation import compute_errors
import argparse
import numpy as np
import os
from tools.coupling_schemes import CouplingScheme


class ProblemType(Enum):
    """
    Enum defines problem type. Details see above.
    """
    DIRICHLET = 1  # Dirichlet problem
    NEUMANN = 2  # Neumann problem


class DomainPart(Enum):
    """
    Enum defines which part of the domain [x_left, x_right] x [y_bottom, y_top] we compute.
    """
    LEFT = 1  # left part of domain
    RIGHT = 2  # right part of domain


class OuterBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and not near(x[0], x_coupling, tol) or near(x[1], y_top, tol) or near(x[1], y_bottom, tol):
            return True
        else:
            return False


class CouplingBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and near(x[0], x_coupling, tol):
            return True
        else:
            return False


def determine_gradient(V_g, u, flux):
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    :param mesh
    :param u: solution where gradient is to be determined
    :return:
    """

    w = TrialFunction(V_g)
    v = TestFunction(V_g)

    a = inner(w, v) * dx
    L = inner(grad(u), v) * dx
    solve(a == L, flux)


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest='dirichlet', action='store_true')
parser.add_argument("-n", "--neumann", help="create a neumann problem", dest='neumann', action='store_true')
parser.add_argument("-wr", "--waveform", nargs=2, default=[1, 1], type=int)
parser.add_argument("-dT", "--window-size", default=1.0, type=float)
parser.add_argument("-cpl", "--coupling-scheme", default=CouplingScheme.SERIAL_FIRST_DIRICHLET.name, type=str)
parser.add_argument("-g", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=0.0, type=float)
parser.add_argument("-tol", "--error-tolerance", help="set accepted error of numerical solution w.r.t analytical solution", default=10**-12, type=float)
parser.add_argument("-dl", "--domain-left", help="right part of the domain is being computed", dest='domain_left', action='store_true')
parser.add_argument("-dr", "--domain-right", help="left part of the domain is being computed", dest='domain_right', action='store_true')
parser.add_argument("-t", "--time-dependence", help="choose whether there is a linear (l), quadratic (q) or sinusoidal (s) dependence on time", type=str, default="l")
parser.add_argument("-mth", "--method", help="time stepping method being used", default='ie')
parser.add_argument("-a", "--arbitrary-coupling-interface", help="uses more general, but less exact method for interpolation on coupling interface, see https://github.com/precice/fenics-adapter/milestone/1", dest='arbitrary_coupling_interface', action='store_true')

args = parser.parse_args()

# coupling parameters
if args.dirichlet:
    problem = ProblemType.DIRICHLET
if args.neumann:
    problem = ProblemType.NEUMANN
if args.dirichlet and args.neumann:
    raise Exception("you can only choose either a dirichlet problem (option -d) or a neumann problem (option -n)")
if not (args.dirichlet or args.neumann):
    raise Exception("you have to choose either a dirichlet problem (option -d) or a neumann problem (option -n)")


# coupling parameters
if args.domain_left:
    domain_part = DomainPart.LEFT
if args.domain_right:
    domain_part = DomainPart.RIGHT
if args.dirichlet and args.neumann:
    raise Exception("you can only choose to either compute the left part of the domain (option -dl) or the right part (option -dr)")
if not (args.domain_left or args.domain_right):
    print("Default domain partitioning is used: Left part of domain is a Dirichlet-type problem; right part is a Neumann-type problem")
    if problem is ProblemType.DIRICHLET:
        domain_part = DomainPart.LEFT
    elif problem is ProblemType.NEUMANN:
        domain_part = DomainPart.RIGHT


# Create mesh and define function space

nx = 10
ny = 10

error_tol = args.error_tolerance

wr_tag = "WR{wr1}{wr2}".format(wr1=args.waveform[0], wr2=args.waveform[1])
window_size = "dT{dT}".format(dT=args.window_size)
coupling_scheme = "{}".format(args.coupling_scheme)
d_subcycling = "D".format(wr_tag=wr_tag)
n_subcycling = "N".format(wr_tag=wr_tag)

configs_path = os.path.join("experiments", wr_tag, window_size, coupling_scheme)

if domain_part is DomainPart.LEFT:
    nx = nx*3
elif domain_part is DomainPart.RIGHT:
    ny = 20

if problem is ProblemType.DIRICHLET:
    adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-D.json")
    other_adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-N.json")
elif problem is ProblemType.NEUMANN:
    adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-N.json")
    other_adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-D.json")

alpha = 3  # parameter alpha
beta = 1.3  # parameter beta
gamma = args.gamma  # parameter gamma, dependence of heat flux on time
y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.5  # x coordinate of coupling interface

if domain_part is DomainPart.LEFT:
    p0 = Point(x_left, y_bottom)
    p1 = Point(x_coupling, y_top)
elif domain_part is DomainPart.RIGHT:
    p0 = Point(x_coupling, y_bottom)
    p1 = Point(x_right, y_top)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 2)

# Define boundary condition
if args.time_dependence == "l":
    print("Linear")
    u_D = Expression('1 + gamma*t*x[0]*x[0] + (1-gamma)*x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
if args.time_dependence == "q":
    print("Quadratic")
    u_D = Expression('1 + gamma * t * t * x[0] * x[0] + (1-gamma) * x[0] * x[0] + alpha * x[1] * x[1] + beta*t', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
elif args.time_dependence == "s":
    print("Sinusoidal")
    u_D = Expression('1 + gamma*sin(t)*x[0]*x[0] + (1-gamma)*x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
u_D_function = interpolate(u_D, V)
# Define flux in x direction on coupling interface (grad(u_D) in normal direction)
if (domain_part is DomainPart.LEFT and problem is ProblemType.DIRICHLET) or \
        (domain_part is DomainPart.RIGHT and problem is ProblemType.NEUMANN):
    if args.time_dependence == "l":
        f_N = Expression('2 * gamma * t * x[0] + 2 * (1-gamma)*x[0] ', degree=1, gamma=gamma, t=0)
    if args.time_dependence == "q":
        f_N = Expression('2 * gamma * t * t * x[0] + 2 * (1-gamma) * x[0] ', degree=1, gamma=gamma, t=0)
    elif args.time_dependence == "s":
        f_N = Expression('2 * gamma * sin(t) * x[0] + 2 * (1-gamma)*x[0] ', degree=1, gamma=gamma, t=0)
elif (domain_part is DomainPart.RIGHT and problem is ProblemType.DIRICHLET) or \
        (domain_part is DomainPart.LEFT and problem is ProblemType.NEUMANN):
    if args.time_dependence == "l":
        f_N = Expression('-2 * gamma * t * x[0] + 2 * (1-gamma)*x[0] ', degree=1, gamma=gamma, t=0)
    if args.time_dependence == "q":
        f_N = Expression('-2 * gamma * t * t * x[0] + 2 * (1-gamma) * x[0] ', degree=1, gamma=gamma, t=0)
    elif args.time_dependence == "s":
        f_N = Expression('-2 * gamma * sin(t) * x[0] + 2 * (1-gamma) * x[0] ', degree=1, gamma=gamma, t=0)
f_N_function = interpolate(f_N, V)

coupling_boundary = CouplingBoundary()
remaining_boundary = OuterBoundary()

bcs = [DirichletBC(V, u_D, remaining_boundary)]
# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

precice = Adapter(adapter_config_filename, other_adapter_config_filename, interpolation_strategy=ExactInterpolationExpression   )  # todo: how to avoid requiring both configs without Waveform Relaxation?

if problem is ProblemType.DIRICHLET:
    dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=u_D_function,
                            write_field=f_N_function, u_n=u_n)
elif problem is ProblemType.NEUMANN:
    dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=f_N_function,
                            write_field=u_D_function, u_n=u_n)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
if args.time_dependence == "l":
    f = Expression('beta + gamma * x[0] * x[0] - 2 * gamma * t - 2 * (1-gamma) - 2 * alpha', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
    f_n = Expression('beta + gamma * x[0] * x[0] - 2 * gamma * t - 2 * (1-gamma) - 2 * alpha', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
if args.time_dependence == "q":
    f = Expression('beta + 2 * gamma * x[0] * x[0] * t - 2 * gamma * t * t- 2 * (1-gamma) - 2 * alpha', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
    f_n = Expression('beta + 2 * gamma * x[0] * x[0] * t - 2 * gamma * t * t - 2 * (1-gamma) - 2 * alpha', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
elif args.time_dependence == "s":
    f = Expression('beta + gamma * x[0] * x[0] * cos(t) - 2 * gamma * sin(t) - 2 * (1-gamma) - 2 * alpha', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
    f_n = Expression('beta + gamma * x[0] * x[0] * cos(t) - 2 * gamma * sin(t) - 2 * (1-gamma) - 2 * alpha', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)

t = 0

if args.method == "ie":
    F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f) * v * dx
elif args.method == "tr":
    F = u * v / dt * dx + dot(grad(u), grad(v)) / 2 * dx + dot(grad(u_n), grad(v)) / 2 * dx - (u_n / dt + f / 2 + f_n / 2) * v * dx

if problem is ProblemType.DIRICHLET:
    # apply Dirichlet boundary condition on coupling interface
    bcs.append(precice.create_coupling_dirichlet_boundary_condition(V, dt))
if problem is ProblemType.NEUMANN:
    # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
    
    if args.method == "tr":
        F += 0.5 * precice.create_coupling_neumann_boundary_condition(v, dt) + 0.5 * precice.create_coupling_neumann_boundary_condition(v, Constant(0))
    elif args.method == "ie":
        F += precice.create_coupling_neumann_boundary_condition(v, dt)

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
u_np1.rename("Temperature", "")

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

temperature_out = File("out/%s.pvd" % precice._solver_name)
ref_out = File("out/ref%s.pvd" % precice._solver_name)
error_out = File("out/error%s.pvd" % precice._solver_name)
flux_out = File("out/flux%s.pvd" % precice._solver_name)

# output solution and reference solution at t=0, n=0
n = 0
temperature_out << u_n
ref_out << u_ref

error_total, error_pointwise = compute_errors(u_n, u_ref, V)
error_out << error_pointwise

# set t_1 = t_0 + dt, this gives u_D^1
u_D.t = t+dt(0)  # call dt(0) to evaluate FEniCS Constant. Todo: is there a better way?
f.t = t+dt(0)
f_n.t = t

V_g = VectorFunctionSpace(mesh, 'P', 1)
flux = Function(V_g)
flux.rename("Flux", "")

determine_gradient(V_g, u_n, flux)
flux_out << flux

while precice.is_coupling_ongoing():

    # Compute solution u^n+1, use bcs u_D^n+1, u^n and coupling bcs
    solve(a == L, u_np1, bcs)

    determine_gradient(V_g, u_np1, flux)

    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
        flux_x, flux_y = flux.split()
        if domain_part is DomainPart.RIGHT:
            flux_x.assign(-flux_x)
        t, n, precice_timestep_complete, precice_dt = precice.advance(flux_x, u_np1, u_n, t, dt(0), n)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem samples temperature on boundary from solution and sends temperature to Dirichlet problem
        t, n, precice_timestep_complete, precice_dt = precice.advance(u_np1, u_np1, u_n, t, dt(0), n)


    dt.assign(np.min([precice.fenics_dt, precice_dt]))  # todo we could also consider deciding on time stepping size inside the adapter

    if precice_timestep_complete:
        u_ref = interpolate(u_D, V)
        u_ref.rename("reference", " ")
        # output solution and reference solution at t_n+1
        error, error_pointwise = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
        temperature_out << u_n
        ref_out << u_ref
        error_out << error_pointwise
        flux_out << flux

    # Update dirichlet BC
    u_D.t = t + dt(0)
    f.t = t + dt(0)
    f_n.t = t

# Hold plot
precice.finalize()
