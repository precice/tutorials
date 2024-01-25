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
from my_enums import ProblemType, DomainPart, TimeDependence
import argparse
import numpy as np
import sympy as sp
from problem_setup import get_geometry, get_manufactured_solution
import pandas as pd
from pathlib import Path


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
parser.add_argument(
    "-s",
    "--n-substeps",
    help="Number of substeps in one window for this participant",
    type=int,
    default=1)
parser.add_argument("-p", "--polynomial-order", help="Polynomial order of manufactured solution", type=int, default=1)
parser.add_argument("-t", "--time-dependence", help="Time dependence of manufactured solution",
                    choices=[e.value for e in TimeDependence], default=TimeDependence.POLYNOMIAL.value)
parser.add_argument("-e", "--error-tol", help="set error tolerance", type=float, default=10**-6,)

args = parser.parse_args()

# Error is bounded by coupling accuracy. In theory we would obtain the analytical solution.
error_tol = args.error_tol

alpha = 3  # parameter alpha
beta = 1.2  # parameter beta

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
u_manufactured, symbols = get_manufactured_solution(TimeDependence(
    args.time_dependence), alpha, beta, p=args.polynomial_order)
u_D = Expression(sp.printing.ccode(u_manufactured), degree=2, t=0)
u_D_function = interpolate(u_D, V)

if problem is ProblemType.DIRICHLET:
    # Define flux in x direction
    flux_x_manufactured = sp.diff(u_manufactured, symbols['x'])
    f_N = Expression(sp.printing.ccode(flux_x_manufactured), degree=1, t=0)
    f_N_function = interpolate(f_N, W)

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

precice, precice_dt, initial_data = None, 0.0, None

# Initialize the adapter according to the specific participant
if problem is ProblemType.DIRICHLET:
    precice = Adapter(adapter_config_filename="precice-adapter-config-D.json")
    precice.initialize(coupling_boundary, read_function_space=V, write_object=f_N_function)
    precice_dt = precice.get_max_time_step_size()
elif problem is ProblemType.NEUMANN:
    precice = Adapter(adapter_config_filename="precice-adapter-config-N.json")
    precice.initialize(coupling_boundary, read_function_space=W, write_object=u_D_function)
    precice_dt = precice.get_max_time_step_size()

dt = Constant(0)
window_dt = precice_dt  # store for later
fenics_dt = precice_dt / args.n_substeps  # use window size provided by preCICE as time step size
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f_manufactured = - sp.diff(sp.diff(u_manufactured,
                                   symbols['x']),
                           symbols['x']) - sp.diff(sp.diff(u_manufactured,
                                                           symbols['y']),
                                                   symbols['y']) + sp.diff(u_manufactured,
                                                                           symbols['t'])
f = Expression(sp.printing.ccode(f_manufactured), degree=2, t=0)
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
    F += v * coupling_expression * ds

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
output_dir = Path("output")
temperature_out = File(str(output_dir / f"{precice.get_participant_name()}.pvd"))
ref_out = File(str(output_dir / f"ref{precice.get_participant_name()}.pvd"))
error_out = File(str(output_dir / f"error{precice.get_participant_name()}.pvd"))
ranks = File(str(output_dir / f"ranks{precice.get_participant_name()}.pvd"))

# output solution and reference solution at t=0, n=0
n = 0
print('output u^%d and u_ref^%d' % (n, n))
temperature_out << u_n
ref_out << u_ref
ranks << mesh_rank

error_total, error_pointwise = compute_errors(u_n, u_ref, V)
error_out << error_pointwise
errors = []
times = []

if problem is ProblemType.DIRICHLET:
    flux = Function(V_g)
    flux.rename("Heat-Flux", "")

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, n)


    precice_dt = precice.get_max_time_step_size()
    dt.assign(np.min([fenics_dt, precice_dt]))

    # Dirichlet BC and RHS need to point to end of current timestep
    u_D.t = t + float(dt)
    f.t = t + float(dt)

    # Coupling BC needs to point to end of current timestep
    read_data = precice.read_data(dt)

    # Update the coupling expression with the new read data
    precice.update_coupling_expression(coupling_expression, read_data)

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

    precice.advance(dt)
    precice_dt = precice.get_max_time_step_size()

    # roll back to checkpoint
    if precice.requires_reading_checkpoint():
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
        errors.append(error)
        times.append(t)

    # Update Dirichlet BC
    u_D.t = t + float(dt)
    f.t = t + float(dt)

precice.finalize()

df = pd.DataFrame()
df["times"] = times
df["errors"] = errors
df = df.set_index('times')
metadata = f'''# time_window_size: {window_dt}
# time_step_size: {fenics_dt}
'''

errors_csv = Path(f"errors-{problem.value}.csv")
errors_csv.unlink(missing_ok=True)

with open(errors_csv, 'a') as f:
    f.write(f"{metadata}")
    df.to_csv(f)
