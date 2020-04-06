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
from fenics import Function, FunctionSpace, Expression, Constant, DirichletBC, TrialFunction, TestFunction, \
    File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate, VectorFunctionSpace
from fenicsadapter import Adapter, InterpolationType
from errorcomputation import compute_errors
from my_enums import ProblemType, Subcycling
import argparse
import numpy as np
from problem_setup import get_geometry, get_problem_setup
import dolfin
from dolfin import FacetNormal, dot

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


parser = argparse.ArgumentParser(description='Solving heat equation for simple or complex interface case')
parser.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest='dirichlet', action='store_true')
parser.add_argument("-n", "--neumann", help="create a neumann problem", dest='neumann', action='store_true')
parser.add_argument("-g", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=0.0, type=float)
parser.add_argument("-a", "--arbitrary-coupling-interface", help="uses more general, but less exact method for interpolation on coupling interface, see https://github.com/precice/fenics-adapter/milestone/1", action='store_true')
parser.add_argument("-i", "--interface", metavar="interface_type string", type=str, choices=['simple', 'complex'], help="Type of coupling interface case to be solved. Options: simple, complex", default="simple")
parser.add_argument("-dom", "--domain", metavar='domain_type string', type=str, choices=['left', 'right', 'circular', 'rectangle'], help="Specifying part of the domain being solved. For simple interface the options are left, right, for complex interface the options are circular, rest")

args = parser.parse_args()

config_file_name = "precice-config.xml"  # TODO should be moved into config, see https://github.com/precice/fenics-adapter/issues/5 , in case file doesnt not exist open will fail

subcycle = Subcycling.NONE

fenics_dt = None
error_tol = None

# for all scenarios, we assume precice_dt == .1
if subcycle is Subcycling.NONE and not args.arbitrary_coupling_interface:
    fenics_dt = .1  # time step size
    error_tol = 10 ** -7  # Error is bounded by coupling accuracy. In theory we would obtain the analytical solution.
    print("Subcycling = NO, Arbitrary coupling interface = NO, error tolerance = {}".format(error_tol))
elif subcycle is Subcycling.NONE and args.arbitrary_coupling_interface:
    fenics_dt = .1  # time step size
    error_tol = 10 ** -3  # error low, if we do not subcycle. In theory we would obtain the analytical solution.
    # TODO For reasons, why we currently still have a relatively high error, see milestone https://github.com/precice/fenics-adapter/milestone/1
    print("Subcycling = NO, Arbitrary coupling interface = YES, error tolerance = {}".format(error_tol))
elif subcycle is Subcycling.MATCHING:
    fenics_dt = .01  # time step size
    error_tol = 10 ** -2  # error increases. If we use subcycling, we cannot assume that we still get the exact solution.
    # TODO Using waveform relaxation, we should be able to obtain the exact solution here, as well.
    print("Subcycling = YES, Matching. error tolerance = {}".format(error_tol))
elif subcycle is Subcycling.NONMATCHING:
    fenics_dt = .03  # time step size
    error_tol = 10 ** -1  # error increases. If we use subcycling, we cannot assume that we still get the exact solution.
    # TODO Using waveform relaxation, we should be able to obtain the exact solution here, as well.
    print("Subcycling = YES, Non-matching. error tolerance = {}".format(error_tol))

alpha = 3  # parameter alpha
beta = 1.3  # parameter beta
gamma = args.gamma  # parameter gamma, dependence of heat flux on time

# Create mesh and separate mesh components for grid, boundary and coupling interface
domain_part, problem = get_problem_setup(args)
mesh, coupling_boundary, remaining_boundary = get_geometry(domain_part)

adapter_config_filename = None
if problem is ProblemType.DIRICHLET:
    adapter_config_filename = "precice-adapter-config-D.json"
elif problem is ProblemType.NEUMANN:
    adapter_config_filename = "precice-adapter-config-N.json"

# Define function space using mesh
V = FunctionSpace(mesh, 'P', 2)
V_g = VectorFunctionSpace(mesh, 'P', 1)

# Define boundary conditions
u_D = Expression('1 + gamma*t*x[0]*x[0] + (1-gamma)*x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
u_D_function = interpolate(u_D, V)
# Define flux in x direction on coupling interface (grad(u_D) in normal direction)
f_N = Expression(("2 * gamma*t*x[0] + 2 * (1-gamma)*x[0]", "2 * alpha*x[1]"), degree=1, gamma=gamma, alpha=alpha, t=0)
f_N_function = interpolate(f_N, V_g)

bcs = [DirichletBC(V, u_D, remaining_boundary)]

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

# Adapter definition and initialization
precice = Adapter(adapter_config_filename)
# Selecting interpolation strategy
precice.set_interpolation_type(InterpolationType.CUBIC_SPLINE)
precice_dt = precice.initialize(coupling_boundary, mesh)

boundary_marker = False
coupling_function = None

# Initialize data to non-standard initial values according to which problem is being solved
if problem is ProblemType.DIRICHLET:
    coupling_function = precice.initialize_data(u_D_function, f_N_function, V)
elif problem is ProblemType.NEUMANN:
    coupling_function = precice.initialize_data(f_N_function, u_D_function, V_g)

# Assigning appropriate dt
dt = Constant(0)
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('beta + gamma * x[0] * x[0] - 2 * gamma * t - 2 * (1-gamma) - 2 * alpha', degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f) * v * dx

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
u_np1.rename("Temperature", "")
t = 0

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

temperature_out = File("out/%s.pvd" % precice.get_solver_name())
ref_out = File("out/ref%s.pvd" % precice.get_solver_name())
error_out = File("out/error%s.pvd" % precice.get_solver_name())

# output solution and reference solution at t=0, n=0
n = 0
print('output u^%d and u_ref^%d' % (n, n))
temperature_out << u_n
ref_out << u_ref

error_total, error_pointwise = compute_errors(u_n, u_ref, V)
error_out << error_pointwise

# set t_1 = t_0 + dt, this gives u_D^1
u_D.t = t + dt(0)  # call dt(0) to evaluate FEniCS Constant. Todo: is there a better way?
f.t = t + dt(0)

flux = Function(V_g)
flux.rename("Flux", "")

while precice.is_coupling_ongoing():

    if precice.is_action_required(precice.action_write_checkpoint()):  # write checkpoint
        precice.store_checkpoint(u_n, t, n)

    # read data from preCICE and get a coupling function
    coupling_function = precice.read()

    # Updating boundary condition at coupling interface
    if problem is ProblemType.DIRICHLET:
        # modify Dirichlet boundary condition on coupling interface
        bcs.append(dolfin.DirichletBC(V, coupling_function, coupling_boundary))
    if problem is ProblemType.NEUMANN:
        # modify Neumann boundary condition on coupling interface, modify weak form correspondingly
        if not boundary_marker:  # there is only 1 Neumann-BC which is at the coupling boundary -> integration over whole boundary
            if coupling_function.is_scalar_valued():
                F += v * coupling_function * dolfin.ds  # this term has to be added to weak form to add a Neumann BC (see e.g. p. 83ff Langtangen, Hans Petter, and Anders Logg. "Solving PDEs in Python The FEniCS Tutorial Volume I." (2016).)
            elif coupling_function.is_vector_valued():
                normal = FacetNormal(mesh)
                F += -v * dot(normal, coupling_function) * dolfin.ds
            else:
                raise Exception("invalid!")
        else:  # For multiple Neumann BCs integration should only be performed over the respective domain.
            # TODO: fix the problem here
            raise Exception("Boundary markers are not implemented yet")
            # return dot(coupling_bc_expression, v) * dolfin.dss(boundary_marker)

    # Assign the correct time step
    dt.assign(np.min([fenics_dt, precice_dt]))

    # Compute solution u^n+1, use bcs u_D^n+1, u^n and coupling bcs
    solve(a == L, u_np1, bcs)

    # Write data to preCICE according to which problem is being solved
    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem reads temperature and writes flux on boundary to Neumann problem
        determine_gradient(V_g, u_np1, flux)
        precice.write(flux.copy())
    elif problem is ProblemType.NEUMANN:
        # Neumann problem reads flux and writes temperature on boundary to Dirichlet problem
        precice.write(u_np1.copy())

    # Call to advance coupling, also returns the optimum time step value
    precice_dt = precice.advance(dt(0))

    # Either revert to old step if timestep has not converged or move to next timestep
    if precice.is_action_required(precice.action_read_checkpoint()):  # roll back to checkpoint
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:
        u_n.assign(u_np1)
        t += dt
        n += 1

    if precice.is_time_window_complete():
        u_ref = interpolate(u_D, V)
        u_ref.rename("reference", " ")
        print('Computing error for: n = %d, t = %.2f' % (n, t))
        error, error_pointwise = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
        print('n = %d, t = %.2f: L2 error on domain = %.3g' % (n, t, error))
        # output solution and reference solution at t_n+1
        print('output u^%d and u_ref^%d' % (n, n))
        temperature_out << u_n
        ref_out << u_ref
        error_out << error_pointwise

    # Update Dirichlet BC
    u_D.t = t + dt(0)
    f.t = t + dt(0)

# Hold plot
precice.finalize()
