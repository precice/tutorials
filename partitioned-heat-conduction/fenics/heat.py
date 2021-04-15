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
    File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate, VectorFunctionSpace, MeshFunction, MPI
from fenicsprecice import Adapter
from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart
import argparse
import numpy as np
from problem_setup import get_geometry
import dolfin
from dolfin import FacetNormal, dot


def determine_gradient(V_g, u, flux):
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    :param V_g: Vector function space
    :param u: solution where gradient is to be determined
    :param flux: returns calculated flux into this value
    """

    w = TrialFunction(V_g)
    v = TestFunction(V_g)

    a = inner(w, v) * dx
    L = inner(grad(u), v) * dx
    solve(a == L, flux)


parser = argparse.ArgumentParser(description="Solving heat equation for simple or complex interface case")
command_group = parser.add_mutually_exclusive_group(required=True)
command_group.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest="dirichlet",
                           action="store_true")
command_group.add_argument("-n", "--neumann", help="create a neumann problem", dest="neumann", action="store_true")
parser.add_argument("-e", "--error-tol", help="set error tolerance", type=float, default=10**-6,)

args = parser.parse_args()

fenics_dt = .1  # time step size
# Error is bounded by coupling accuracy. In theory we would obtain the analytical solution.
error_tol = args.error_tol

alpha = 3  # parameter alpha
beta = 1.3  # parameter beta

if args.dirichlet and not args.neumann:
    problem = ProblemType.DIRICHLET
    domain_part = DomainPart.LEFT
elif args.neumann and not args.dirichlet:
    problem = ProblemType.NEUMANN
    domain_part = DomainPart.RIGHT

mesh, coupling_boundary, remaining_boundary = get_geometry(domain_part)

# Define function space using mesh
V = FunctionSpace(mesh, 'P', 2)
V_g = VectorFunctionSpace(mesh, 'P', 1)
W = V_g.sub(0).collapse()

# Define boundary conditions
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0)
u_D_function = interpolate(u_D, V)

if problem is ProblemType.DIRICHLET:
    # Define flux in x direction
    f_N = Expression("2 * x[0]", degree=1, alpha=alpha, t=0)
    f_N_function = interpolate(f_N, W)

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

precice, precice_dt, initial_data = None, 0.0, None

# Initialize the adapter according to the specific participant
if problem is ProblemType.DIRICHLET:
    precice = Adapter(adapter_config_filename="precice-adapter-config-D.json")
    precice_dt = precice.initialize(coupling_boundary, read_function_space=V, write_object=f_N_function)
elif problem is ProblemType.NEUMANN:
    precice = Adapter(adapter_config_filename="precice-adapter-config-N.json")
    precice_dt = precice.initialize(coupling_boundary, read_function_space=W, write_object=u_D_function)

dt = Constant(0)
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('beta - 2 - 2*alpha', degree=2, alpha=alpha, beta=beta, t=0)
F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f) * v * dx

bcs = [DirichletBC(V, u_D, remaining_boundary)]

# Set boundary conditions at coupling interface once wrt to the coupling
# expression
coupling_expression = precice.create_coupling_expression()
if problem is ProblemType.DIRICHLET:
    # modify Dirichlet boundary condition on coupling interface
    bcs.append(DirichletBC(V, coupling_expression, coupling_boundary))
if problem is ProblemType.NEUMANN:
    # modify Neumann boundary condition on coupling interface, modify weak
    # form correspondingly
    F += v * coupling_expression * dolfin.ds

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
u_np1.rename("Temperature", "")
t = 0

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

# mark mesh w.r.t ranks
mesh_rank = MeshFunction("size_t", mesh, mesh.topology().dim())
if problem is ProblemType.NEUMANN:
    mesh_rank.set_all(MPI.rank(MPI.comm_world) + 4)
else:
    mesh_rank.set_all(MPI.rank(MPI.comm_world) + 0)
mesh_rank.rename("myRank", "")

# Generating output files
temperature_out = File("out/%s.pvd" % precice.get_participant_name())
ref_out = File("out/ref%s.pvd" % precice.get_participant_name())
error_out = File("out/error%s.pvd" % precice.get_participant_name())
ranks = File("out/ranks%s.pvd" % precice.get_participant_name())

# output solution and reference solution at t=0, n=0
n = 0
print('output u^%d and u_ref^%d' % (n, n))
temperature_out << u_n
ref_out << u_ref
ranks << mesh_rank

error_total, error_pointwise = compute_errors(u_n, u_ref, V)
error_out << error_pointwise

# set t_1 = t_0 + dt, this gives u_D^1
# call dt(0) to evaluate FEniCS Constant. Todo: is there a better way?
u_D.t = t + dt(0)
f.t = t + dt(0)

if problem is ProblemType.DIRICHLET:
    flux = Function(V_g)
    flux.rename("Flux", "")

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.is_action_required(precice.action_write_iteration_checkpoint()):
        precice.store_checkpoint(u_n, t, n)

    read_data = precice.read_data()

    # Update the coupling expression with the new read data
    precice.update_coupling_expression(coupling_expression, read_data)

    dt.assign(np.min([fenics_dt, precice_dt]))

    # Compute solution u^n+1, use bcs u_D^n+1, u^n and coupling bcs
    solve(a == L, u_np1, bcs)

    # Write data to preCICE according to which problem is being solved
    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem reads temperature and writes flux on boundary to Neumann problem
        determine_gradient(V_g, u_np1, flux)
        flux_x = interpolate(flux.sub(0), W)
        precice.write_data(flux_x)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem reads flux and writes temperature on boundary to Dirichlet problem
        precice.write_data(u_np1)

    precice_dt = precice.advance(dt(0))

    # roll back to checkpoint
    if precice.is_action_required(precice.action_read_iteration_checkpoint()):
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:  # update solution
        u_n.assign(u_np1)
        t += float(dt)
        n += 1

    if precice.is_time_window_complete():
        u_ref = interpolate(u_D, V)
        u_ref.rename("reference", " ")
        error, error_pointwise = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
        print('n = %d, t = %.2f: L2 error on domain = %.3g' % (n, t, error))
        # output solution and reference solution at t_n+1
        print('output u^%d and u_ref^%d' % (n, n))
        temperature_out << u_n
        ref_out << u_ref
        error_out << error_pointwise

    # Update Dirichlet BC
    u_D.t = t + float(dt)
    f.t = t + float(dt)

# Hold plot
precice.finalize()
