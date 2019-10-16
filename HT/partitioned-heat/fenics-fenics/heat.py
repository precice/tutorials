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
from fenics import Function, FunctionSpace, Constant, DirichletBC, \
    TrialFunction, TestFunction, File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate, VectorFunctionSpace, set_log_level
from fenicsadapter import Adapter, ExactInterpolationExpression, GeneralInterpolationExpression
from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart
from problem_setup import get_problem_setup, get_geometry, OuterBoundary, CouplingBoundary, CompleteBoundary
import argparse
import numpy as np
import os
import sdc.simple_sdc
from tools.coupling_schemes import CouplingScheme

set_log_level(30)

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
parser.add_argument("-wri", "--waveform-interpolation-strategy", help="specify interpolation strategy used by waveform relaxation", default="linear", choices=['linear', 'quadratic', 'cubic', 'quartic'], type=str)
parser.add_argument("-dT", "--window-size", default=1.0, type=float)
parser.add_argument("-T", "--total-time", default=10, type=float)
parser.add_argument("-cpl", "--coupling-scheme", default=CouplingScheme.SERIAL_FIRST_DIRICHLET.name, type=str)
parser.add_argument("-gamma", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=1.0, type=float)
parser.add_argument("-alpha", "--alpha", help="parameter gamma to set temporal dependence of heat flux", default=3.0, type=float)
parser.add_argument("-beta", "--beta", help="parameter gamma to set temporal dependence of heat flux", default=1.2, type=float)
parser.add_argument("-tol", "--error-tolerance", help="set accepted error of numerical solution w.r.t analytical solution", default=10**-12, type=float)
parser.add_argument("-dl", "--domain-left", help="right part of the domain is being computed", dest='domain_left', action='store_true')
parser.add_argument("-dr", "--domain-right", help="left part of the domain is being computed", dest='domain_right', action='store_true')
parser.add_argument("-t", "--time-dependence", help="choose whether there is a linear (l), quadratic (q), cubic (c) or sinusoidal (s) dependence on time", type=str, default="l", choices=['l', 'q', 'c', 'qrt', 's'])
parser.add_argument("-mth", "--method", help="time stepping method being used", default='ie', choices=['ie', 'tr', 'sdc'])
parser.add_argument("--sdc-K", help="number of correction sweeps used for SDC", type=int, default=16)
parser.add_argument("-nx", "--nx", help="number of DoFs in x direction", type=int, default=20)
parser.add_argument("-ny", "--ny", help="number of DoFs in y direction", type=int, default=20)
parser.add_argument("-a", "--arbitrary-coupling-interface", help="uses more general, but less exact method for interpolation on coupling interface, see https://github.com/precice/fenics-adapter/milestone/1", dest='arbitrary_coupling_interface', action='store_true')
parser.add_argument("-ovtk", "--output-vtk", help="provide vtk output", action='store_true')
parser.add_argument("-m", "--monolithic", help="switch to monolithic case", action='store_true')

args = parser.parse_args()

if not args.monolithic:
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
else:
    domain_part = DomainPart.ALL
    problem = ProblemType.DIRICHLET

error_tol = args.error_tolerance

wr_tag = "WR{wr1}{wr2}".format(wr1=args.waveform[0], wr2=args.waveform[1])
window_size = "dT{dT}".format(dT=args.window_size)
coupling_scheme = "{}".format(args.coupling_scheme)
d_subcycling = "D".format(wr_tag=wr_tag)
n_subcycling = "N".format(wr_tag=wr_tag)

t_total = args.total_time

if not args.monolithic:
    configs_path = os.path.join("experiments", wr_tag, window_size, coupling_scheme)

    if problem is ProblemType.DIRICHLET:
        adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-D.json")
        other_adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-N.json")
    elif problem is ProblemType.NEUMANN:
        adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-N.json")
        other_adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-D.json")

# Create mesh and define function space
nx = args.nx
ny = args.ny
if args.monolithic:
    nx *= 2

mesh = get_geometry(domain_part, nx, ny)
V = FunctionSpace(mesh, 'P', 2)

# Get Expressions defining boundary conditions, right hand side and analytical solution of the problem
f_np1, f_n, u_D, f_N = get_problem_setup(args, domain_part, problem)

u_D_function = interpolate(u_D, V)
f_N_function = interpolate(f_N, V)

if args.monolithic:
    complete_boundary = CompleteBoundary()
    bcs = [DirichletBC(V, u_D, complete_boundary)]
else:
    coupling_boundary = CouplingBoundary()
    remaining_boundary = OuterBoundary()
    bcs = [DirichletBC(V, u_D, remaining_boundary)]

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

if not args.monolithic:
    precice = Adapter(adapter_config_filename, other_adapter_config_filename, interpolation_strategy=ExactInterpolationExpression, wr_interpolation_strategy=args.waveform_interpolation_strategy)  # todo: how to avoid requiring both configs without Waveform Relaxation?

    if problem is ProblemType.DIRICHLET:
        dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=u_D_function, write_field=f_N_function, u_n=u_n)
    elif problem is ProblemType.NEUMANN:
        dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=f_N_function, write_field=u_D_function, u_n=u_n)
else:
    dt = Constant(args.window_size)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

t = 0

if args.method == "ie":
    F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f_np1) * v * dx
elif args.method == "tr":
    F = u * v / dt * dx + dot(grad(u), grad(v)) / 2 * dx + dot(grad(u_n), grad(v)) / 2 * dx - (u_n / dt + f_np1 / 2 + f_n / 2) * v * dx
elif args.method == "sdc":
    times = sdc.simple_sdc.x * dt
    dts = times[1:] - times[:-1]
    F = []
    for i in range(dts.size):
        F_i = u * v / dts[i] * dx + dot(grad(u), grad(v)) * dx - (u_n / dts[i] + f_np1) * v * dx
        F.append(F_i)

if args.method == "sdc":
    bc_non_coupling = bcs[0]
    bcs = []
    for i in range(times.size):
        bcs.append([])
        bcs[i].append(bc_non_coupling)

if not args.monolithic:
    if problem is ProblemType.DIRICHLET:
        # apply Dirichlet boundary condition on coupling interface
        if args.method == "tr" or args.method == "ie":
            bcs.append(precice.create_coupling_dirichlet_boundary_condition(V, dt))
        elif args.method == "sdc":
            for i in range(times.size):
                bcs[i].append(precice.create_coupling_dirichlet_boundary_condition(V, times[i]))

    if problem is ProblemType.NEUMANN:
        # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
    
        if args.method == "tr":
            F += 0.5 * precice.create_coupling_neumann_boundary_condition(v, dt) + 0.5 * precice.create_coupling_neumann_boundary_condition(v, Constant(0))
        elif args.method == "ie":
            F += precice.create_coupling_neumann_boundary_condition(v, dt)
        elif args.method == "sdc":
            for i in range(dts.size):
                F[i] += precice.create_coupling_neumann_boundary_condition(v, times[i + 1])

# Time-stepping
u_np1 = Function(V)
u_np1.rename("Temperature", "")

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

# output solution and reference solution at t=0, n=0
n = 0

error_total, error_pointwise = compute_errors(u_n, u_ref, V)


# set t_1 = t_0 + dt, this gives u_D^1
u_D.t = t+dt(0)  # call dt(0) to evaluate FEniCS Constant. Todo: is there a better way?
f_np1.t = t + dt(0)
f_n.t = t

V_g = VectorFunctionSpace(mesh, 'P', 1)
flux = Function(V_g)
flux.rename("Flux", "")

determine_gradient(V_g, u_n, flux)

if args.output_vtk:
    temperature_out = File("out/%s.pvd" % problem.name)
    ref_out = File("out/ref%s.pvd" % problem.name)
    error_out = File("out/error%s.pvd" % problem.name)
    flux_out = File("out/flux%s.pvd" % problem.name)
    temperature_out << u_n
    ref_out << u_ref
    error_out << error_pointwise
    flux_out << flux

# only needed for SDC


def implicit_euler(y0, t0, dt_step, weak_form, bc_expr, bc):
    u_n.assign(y0)
    bc_expr.t = t0 + dt_step
    f_np1.t = t0 + dt_step
    y1 = Function(V)
    solve(lhs(weak_form) == rhs(weak_form), y1, bc)
    return y1


def black_box_implicit_euler(y0, t0, dt_step, i):
    r_val = implicit_euler(y0, t0, dt_step, F[i], u_D, bcs[i+1])
    return r_val


def compute_rhs(y0, t0, i):
    from fenics import div
    w = TrialFunction(V)
    v = TestFunction(V)
    u_D.t = t0
    f_n.t = t0  # does this interfere with use of f at other places?
    u_D_function = interpolate(u_D, V)
    a = w * v * dx
    L = (div(grad(y0)) + f_n) * v * dx
    u_rhs = Function(V)
    solve(a == L, u_rhs, bcs[i])
    return u_rhs

# ending here

time_loop = True
precice_timestep_complete = False

if args.monolithic:
    time_loop = t < t_total
else:
    time_loop = precice.is_coupling_ongoing()

while time_loop:

    # Compute solution u^n+1, use bcs u_D^n+1, u^n and coupling bcs
    if args.method == 'ie' or args.method =='tr':
        solve(lhs(F) == rhs(F), u_np1, bcs)
    elif args.method == 'sdc':
        u_np1 = sdc.simple_sdc.sdc_step(u_n, t, black_box_implicit_euler, compute_rhs, dt(0), V, K=args.sdc_K)

    if not args.monolithic:
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
    else:
        u_n.assign(u_np1)
        t += dt(0)
        n += 1

    if precice_timestep_complete or args.monolithic:
        u_ref = interpolate(u_D, V)
        u_ref.rename("reference", " ")
        # output solution and reference solution at t_n+1
        error, error_pointwise = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
        if args.output_vtk:
            temperature_out << u_n
            ref_out << u_ref
            error_out << error_pointwise
            flux_out << flux

    # Update dirichlet BC
    u_D.t = t + dt(0)
    f_np1.t = t + dt(0)
    f_n.t = t

    if args.monolithic:
        time_loop = t < t_total
    else:
        time_loop = precice.is_coupling_ongoing()

print(t)
u_D.t = t
u_ref = interpolate(u_D, V)
error_last, _ = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
print(error_last)
print(dt(0))

if not args.monolithic:
    precice.finalize()

with open("errors_{scheme}_{total_time}_{participant}".format(scheme=args.method, total_time=t, participant=problem.name), 'a') as outfile:
    outfile.write("{}, {}\n".format(dt(0), error_last))
