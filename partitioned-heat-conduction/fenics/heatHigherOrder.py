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

For information on the partitioned heat conduction problem using higher-order implicit Runge-Kutta methods see
* "Wullenweber, Nikola. High-order time stepping schemes for solving partial differential equations with FEniCS. Bachelor's thesis at Technical University of Munich, 2021. URL: https://mediatum.ub.tum.de/1621360"
* "Vinnitchenko, Niklas. Evaluation of higher-order coupling schemes with FEniCS-preCICE. Bachelor's thesis at Technical University of Munich, 2023. URL: https://mediatum.ub.tum.de/1732367"

The implementation of the higher-order implicit Runge-Kutta methods is based on: https://github.com/NikoWul/FenicsIrksome
"""

from __future__ import print_function, division
from fenics import Function, FunctionSpace, Expression, Constant, DirichletBC, TrialFunction, TestFunction, \
    File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate, VectorFunctionSpace, MeshFunction, MPI, MixedElement, split, project
from fenicsprecice import Adapter
from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart, TimeDependence
import argparse
import numpy as np
from problem_setup import get_geometry, get_manufactured_solution
import sympy as sp
from utils.ButcherTableaux import *
import utils.utils as utl
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

# time stepping setup
# scheme
tsm = GaussLegendre(2)
# depending on tsm, we define the trial and test function space
if tsm.num_stages == 1:
    Vbig = V
else:
    # for multi-stage RK methods, we need more dimensional function spaces
    mixed = MixedElement(tsm.num_stages * [V.ufl_element()])
    Vbig = FunctionSpace(V.mesh(), mixed)

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

# stage times of the time stepping scheme relative to the current time and dependent on the current dt
stage_times = [tsm.c[i] * dt for i in range(tsm.num_stages)]

# Define variational problem
# trial and test functions
u = TrialFunction(Vbig)
v = TestFunction(Vbig)
# if dim(Vbig)>1, f needs to be stored in an array with different times,
# because in each stage, of an RK method it is evaluated at a different time
# du_dt-Laplace(u) = f
f_manufactured = - sp.diff(sp.diff(u_manufactured,
                                   symbols['x']),
                           symbols['x']) - sp.diff(sp.diff(u_manufactured,
                                                           symbols['y']),
                                                   symbols['y']) + sp.diff(u_manufactured,
                                                                           symbols['t'])
# initial time assumed to be 0
f = [Expression(sp.ccode(f_manufactured), degree=2, t=float(stage_times[i])) for i in range(tsm.num_stages)]
# get variational form of the problem
F = utl.get_variational_problem(v=v, initial_condition=u_n, dt=dt, fs=f, tsm=tsm, k=u)

# boundary conditions

# we define variational form for each RK stage. According to
# https://doi.org/10.1145/3466168, All of those stages
# represent time derivatives and not the actual solution.
# Thus, we need time derivatives of the Dirichlet boundaries as boundary conditions

# get time derivative of u
du_dt_expr = u_manufactured.diff(symbols["t"])
du_dt = [Expression(sp.ccode(du_dt_expr), degree=2, t=float(stage_times[i])) for i in range(tsm.num_stages)]
# set up boundary conditions
bc = []
# Create for each dimension of Vbig a coupling expression for either the time derivatives (Dirichlet side)
#  or Neumann side (no time derivatives required as they are enforced with changing F)
coupling_expressions = [precice.create_coupling_expression() for _ in range(tsm.num_stages)]

# wrap Vbig into array independent of the length to allow equal treatment below
if tsm.num_stages == 1:
    Vbigs = [Vbig]
    vs = [v]
else:
    Vbigs = [Vbig.sub(i) for i in range(tsm.num_stages)]
    vs = split(v)

# for the boundary which is not the coupling boundary, we can just use the boundary conditions as usual.
# Each stage needs a boundary condition
for i in range(tsm.num_stages):
    if problem is ProblemType.DIRICHLET:
        bc.append(DirichletBC(Vbigs[i], du_dt[i], remaining_boundary))
        bc.append(DirichletBC(Vbigs[i], coupling_expressions[i], coupling_boundary))
    else:
        bc.append(DirichletBC(Vbigs[i], du_dt[i], remaining_boundary))
        F += vs[i] * coupling_expressions[i] * ds

a, L = lhs(F), rhs(F)

# we additionally need the values of the stages from the RK method.
# they are stored here
k = Function(Vbig)

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
print("output u^%d and u_ref^%d" % (n, n))
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
    # update boundary conditions
    for i in range(tsm.num_stages):
        f[i].t = t + float(stage_times[i])
        du_dt[i].t = t + float(stage_times[i])

    # can be deleted, if preCICE offers an API call for getting the time derivative of the waveform
    # (see https://github.com/precice/precice/issues/1908)
    # compute time derivatives for dirichlet side
    if problem is ProblemType.DIRICHLET:
        # approximate the function which preCICE uses with BSplines
        bsplns = utl.b_splines(precice, 5, float(dt))

        # get first derivative
        bsplns_der = {}
        for ki in bsplns.keys():
            bsplns_der[ki] = bsplns[ki].derivative(1)

    if problem is ProblemType.DIRICHLET:
        # if preCICE offers API to get the time derivative directly from it, you can remove the entire for-loop
        # and directly insert it into update_coupling_expression
        for i in range(tsm.num_stages):
            # values of derivative at current time
            val = {}
            for ki in bsplns_der.keys():
                val[ki] = bsplns_der[ki](float(stage_times[i]))
            precice.update_coupling_expression(coupling_expressions[i], val)
    else:
        # Neumann boundaries just need temperature flux
        for i in range(tsm.num_stages):
            precice.update_coupling_expression(coupling_expressions[i], precice.read_data(stage_times[i]))

    # getting the solution of the current time step

    # instead of directly solving for u, we look for the values of k.
    # with those we can assemble the solution and thus get the discrete evolution
    solve(a == L, k, bc)

    if tsm.num_stages == 1:
        ks = [k]
    else:
        ks = [k.sub(i) for i in range(tsm.num_stages)]

    # now we need to add up the stages k according to the time stepping scheme
    # -> assembly of discrete evolution
    u_np1 = u_n
    for i in range(tsm.num_stages):
        u_np1 += dt * tsm.b[i] * ks[i]

    # u is in function space V and not Vbig -> project u_np1 to V
    u_np1 = project(u_np1, V)

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
        print("n = %d, t = %.2f: L2 error on domain = %.3g" % (n, t, error))
        # output solution and reference solution at t_n+1
        print('output u^%d and u_ref^%d' % (n, n))
        temperature_out << u_n
        ref_out << u_ref
        error_out << error_pointwise
        errors.append(error)
        times.append(t)

    # Update Dirichlet BC
    u_D.t = t + float(dt)
    for i in range(tsm.num_stages):
        f[i].t = t + float(stage_times[i])
        du_dt[i].t = t + float(stage_times[i])

# Hold plot
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