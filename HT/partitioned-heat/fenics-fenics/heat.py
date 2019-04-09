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
    TrialFunction, TestFunction, File, solve, plot, lhs, rhs, grad, inner, dot, dx, ds, assemble, interpolate, project, near, VectorFunctionSpace
from enum import Enum
from fenicsadapter import Adapter
from errorcomputation import compute_errors
import argparse
import numpy as np


class ProblemType(Enum):
    """
    Enum defines problem type. Details see above.
    """
    DIRICHLET = 1  # Dirichlet problem
    NEUMANN = 2  # Neumann problem

class Subcycling(Enum):
    """
    Enum defines which kind of subcycling is used
    """
    NONE = 0  # no subcycling, precice_dt == fenics_dt
    MATCHING = 1  # subcycling, where fenics_dt fits into precice_dt, mod(precice_dt, fenics_dt) == 0
    NONMATCHING = 2  # subcycling, where fenics_dt does not fit into precice_dt, mod(precice_dt, fenics_dt) != 0
    DIFFERENT = 3  # subcycling, where fenics_dt fits into precice_dt and fenics_dt differs for the two subdomains

    # note: the modulo expressions above should be understood in an exact way (no floating point round off problems. For
    # details, see https://stackoverflow.com/questions/14763722/python-modulo-on-floats)


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


def fluxes_from_temperature_full_domain(u, V, mesh):
    """
    compute flux from weak form (see p.3 in Toselli, Andrea, and Olof Widlund. Domain decomposition methods-algorithms and theory. Vol. 34. Springer Science & Business Media, 2006.)
    :param u: known temperature field
    :param V: function space
    :param mesh: the underlying mesh
    :return:
    """
    degree = V.ufl_element().degree()
    W = VectorFunctionSpace(mesh, 'P', degree)
    grad_u_x, grad_u_y = project(grad(u), W).split()
    return grad_u_x  # todo: this is not general! In the end we will need the normal flux.


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest='dirichlet', action='store_true')
parser.add_argument("-n", "--neumann", help="create a neumann problem", dest='neumann', action='store_true')
parser.add_argument("-wr", "--waveform", action='store_true')

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

# Create mesh and define function space

nx = 10
ny = 10

if args.waveform:
    subcycle = Subcycling.DIFFERENT
else:
    subcycle = Subcycling.NONE

# for all scenarios, we assume precice_dt == .1
if subcycle is Subcycling.NONE:
    fenics_dt = .1  # time step size
    wr_tag = "WR11"
    d_subcycling = "D-{wr_tag}".format(wr_tag=wr_tag)
    n_subcycling = "N-{wr_tag}".format(wr_tag=wr_tag)
    error_tol = 10 ** -4  # error low, if we do not subcycle. In theory we would obtain the analytical solution.
elif subcycle is Subcycling.MATCHING:
    fenics_dt = .1  # time step size
    wr_tag = "WR22"
    d_subcycling = "D-{wr_tag}".format(wr_tag=wr_tag)
    n_subcycling = "N-{wr_tag}".format(wr_tag=wr_tag)
    error_tol = 10 ** -2  # error increases. If we use subcycling, we cannot assume that we still get the exact solution.
    # TODO Using waveform relaxation, we should be able to obtain the exact solution here, as well.
elif subcycle is Subcycling.NONMATCHING:
    fenics_dt = .03  # time step size
    error_tol = 10 ** -1  # error increases. If we use subcycling, we cannot assume that we still get the exact solution.
    # TODO Using waveform relaxation, we should be able to obtain the exact solution here, as well.
elif subcycle is Subcycling.DIFFERENT:
    if problem is ProblemType.DIRICHLET:
        fenics_dt = .1  # time step size
    elif problem is ProblemType.NEUMANN:
        fenics_dt = .05  # time step size
    error_tol = 10 ** -4
    wr_tag = "WR12"
    d_subcycling = "D-{wr_tag}".format(wr_tag=wr_tag)
    n_subcycling = "N-{wr_tag}".format(wr_tag=wr_tag)
    # TODO Using waveform relaxation, we should be able to obtain the exact solution here, as well.

if problem is ProblemType.DIRICHLET:
    nx = nx*3
    adapter_config_filename = "precice-adapter-config-{d_subcycling}.json".format(d_subcycling=d_subcycling)
    other_adapter_config_filename = "precice-adapter-config-{n_subcycling}.json".format(n_subcycling=n_subcycling)

elif problem is ProblemType.NEUMANN:
    adapter_config_filename = "precice-adapter-config-{n_subcycling}.json".format(n_subcycling=n_subcycling)
    other_adapter_config_filename = "precice-adapter-config-{d_subcycling}.json".format(d_subcycling=d_subcycling)

alpha = 3  # parameter alpha
beta = 1.3  # parameter beta
y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.5  # x coordinate of coupling interface

if problem is ProblemType.DIRICHLET:
    p0 = Point(x_left, y_bottom)
    p1 = Point(x_coupling, y_top)
elif problem is ProblemType.NEUMANN:
    p0 = Point(x_coupling, y_bottom)
    p1 = Point(x_right, y_top)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 2)

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0)
u_D_function = interpolate(u_D, V)
# Define flux in x direction on coupling interface (grad(u_D) in normal direction)
f_N_function = fluxes_from_temperature_full_domain(u_D_function, V, mesh)

coupling_boundary = CouplingBoundary()
remaining_boundary = ComplementaryBoundary(coupling_boundary)

bcs = [DirichletBC(V, u_D, remaining_boundary)]
# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

precice = Adapter(adapter_config_filename, other_adapter_config_filename)  # todo: how to avoid requiring both configs without Waveform Relaxation?

if problem is ProblemType.DIRICHLET:
    precice_dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=u_D_function,
                                    write_field=f_N_function, u_n=u_n)
elif problem is ProblemType.NEUMANN:
    precice_dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=f_N_function,
                                    write_field=u_D_function, u_n=u_n)

dt = Constant(0)
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2 * alpha)
F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f) * v * dx

if problem is ProblemType.DIRICHLET:
    # apply Dirichlet boundary condition on coupling interface
    bcs.append(precice.create_coupling_dirichlet_boundary_condition(V))
if problem is ProblemType.NEUMANN:
    # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
    F += precice.create_coupling_neumann_boundary_condition(v)

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
u_np1.rename("Temperature", "")
t = 0

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

temperature_out = File("out/%s.pvd" % precice._solver_name)
ref_out = File("out/ref%s.pvd" % precice._solver_name)
error_out = File("out/error%s.pvd" % precice._solver_name)

# output solution and reference solution at t=0, n=0
n = 0
print('output u^%d and u_ref^%d' % (n, n))
temperature_out << u_n
ref_out << u_ref

error_total, error_pointwise = compute_errors(u_n, u_ref, V)
error_out << error_pointwise

# set t_1 = t_0 + dt, this gives u_D^1
u_D.t = t + dt(0)  # call dt(0) to evaluate FEniCS Constant. Todo: is there a better way?

while precice.is_coupling_ongoing():

    x_check, y_check = 1.5, 0.5

    if problem is ProblemType.DIRICHLET:
        u_ref = interpolate(u_D, V)
        print("before solve:")
        print(u_n(x_check, y_check))
        print(u_np1(x_check, y_check))
        print(u_ref(x_check, y_check))

    # Compute solution u^n+1, use bcs u_D^n+1, u^n and coupling bcs
    solve(a == L, u_np1, bcs)

    print("t={t}; dt={dt}".format(t=t, dt=dt(0)))

    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
        fluxes = fluxes_from_temperature_full_domain(u_np1, V, mesh)
        t, n, precice_timestep_complete, precice_dt = precice.advance(fluxes, u_np1, u_n, t, dt(0), n)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem samples temperature on boundary from solution and sends temperature to Dirichlet problem
        t, n, precice_timestep_complete, precice_dt = precice.advance(u_np1, u_np1, u_n, t, dt(0), n)

    if problem is ProblemType.DIRICHLET:
        print("after solve:")
        print(u_n(x_check, y_check))
        print(u_np1(x_check, y_check))
        print(u_ref(x_check, y_check))
        print(fluxes(x_check, y_check))

    dt.assign(np.min([fenics_dt, precice_dt]))  # todo we could also consider deciding on time stepping size inside the adapter

    if precice_timestep_complete:
        u_ref = interpolate(u_D, V)
        u_ref.rename("reference", " ")
        error, error_pointwise = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
        print('n = %d, t = %.2f: L2 error on domain = %.3g' % (n, t, error))
        # output solution and reference solution at t_n+1
        print('output u^%d and u_ref^%d' % (n, n))
        temperature_out << u_n
        ref_out << u_ref
        error_out << error_pointwise

    # Update dirichlet BC
    u_D.t = t + dt(0)

# Hold plot
precice.finalize()
