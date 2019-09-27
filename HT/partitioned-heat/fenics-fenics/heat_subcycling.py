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
    TrialFunction, TestFunction, File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate, VectorFunctionSpace
from fenicsadapter import Adapter, ExactInterpolationExpression, GeneralInterpolationExpression
from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart
from problem_setup import get_problem_setup, get_geometry, OuterBoundary, CouplingBoundary
import argparse
import numpy as np
import os
from tools.coupling_schemes import CouplingScheme


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
parser.add_argument("-wri", "--waveform-interpolation-strategy", help="specify interpolation strategy used by waveform relaxation", default="linear", choices=['linear', 'quadratic', 'cubic'], type=str)  # this is a don't-care argument
parser.add_argument("-dT", "--window-size", default=1.0, type=float)
parser.add_argument("-cpl", "--coupling-scheme", default=CouplingScheme.SERIAL_FIRST_DIRICHLET.name, type=str)
parser.add_argument("-gamma", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=1.0, type=float)
parser.add_argument("-alpha", "--alpha", help="parameter gamma to set temporal dependence of heat flux", default=3.0, type=float)
parser.add_argument("-beta", "--beta", help="parameter gamma to set temporal dependence of heat flux", default=0.0, type=float)
parser.add_argument("-tol", "--error-tolerance", help="set accepted error of numerical solution w.r.t analytical solution", default=10**12, type=float)
parser.add_argument("-dl", "--domain-left", help="right part of the domain is being computed", dest='domain_left', action='store_true')
parser.add_argument("-dr", "--domain-right", help="left part of the domain is being computed", dest='domain_right', action='store_true')
parser.add_argument("-t", "--time-dependence", help="choose whether there is a linear (l), quadratic (q) or sinusoidal (s) dependence on time", type=str, default="s")
parser.add_argument("-mth", "--method", help="time stepping method being used", default='ie')
parser.add_argument("-nx", "--nx", help="number of DoFs in x direction", type=int, default=20)
parser.add_argument("-ny", "--ny", help="number of DoFs in y direction", type=int, default=20)
parser.add_argument("-a", "--arbitrary-coupling-interface", help="uses more general, but less exact method for interpolation on coupling interface, see https://github.com/precice/fenics-adapter/milestone/1", dest='arbitrary_coupling_interface', action='store_true')
parser.add_argument("-ovtk", "--output-vtk", help="provide vtk output", action='store_true')

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


error_tol = args.error_tolerance

wr_tag = "WR{wr1}{wr2}".format(wr1=args.waveform[0], wr2=args.waveform[1])
window_size = "dT{dT}".format(dT=args.window_size)
coupling_scheme = "{}".format(args.coupling_scheme)
d_subcycling = "D".format(wr_tag=wr_tag)
n_subcycling = "N".format(wr_tag=wr_tag)

configs_path = os.path.join("experiments", wr_tag, window_size, coupling_scheme)

if problem is ProblemType.DIRICHLET:
    adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-D.json")
    other_adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-N.json")
elif problem is ProblemType.NEUMANN:
    adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-N.json")
    other_adapter_config_filename = os.path.join(configs_path, "precice-adapter-config-D.json")

# Create mesh and define function space
mesh = get_geometry(domain_part, args.nx, args.ny)
V = FunctionSpace(mesh, 'P', 2)

# Get Expressions defining boundary conditions, right hand side and analytical solution of the problem
f_np1, f_n, u_D, f_N = get_problem_setup(args, domain_part, problem)

u_D_function = interpolate(u_D, V)
f_N_function = interpolate(f_N, V)

coupling_boundary = CouplingBoundary()
remaining_boundary = OuterBoundary()

bcs = [DirichletBC(V, u_D, remaining_boundary)]
# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

precice = Adapter(adapter_config_filename, interpolation_strategy=ExactInterpolationExpression   )  # todo: how to avoid requiring both configs without Waveform Relaxation?

if problem is ProblemType.DIRICHLET:
    precice_dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=u_D_function,
                            write_field=f_N_function, u_n=u_n)
elif problem is ProblemType.NEUMANN:
    precice_dt = precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=f_N_function,
                            write_field=u_D_function, u_n=u_n)

if domain_part is DomainPart.LEFT:
    fenics_dt = precice_dt / args.waveform[0]
elif domain_part is DomainPart.RIGHT:
    fenics_dt = precice_dt / args.waveform[1]

dt = Constant(np.min([fenics_dt, precice_dt]))
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

t = 0

if args.method == "ie":
    F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f_np1) * v * dx
elif args.method == "tr":
    F = u * v / dt * dx + dot(grad(u), grad(v)) / 2 * dx + dot(grad(u_n), grad(v)) / 2 * dx - (u_n / dt + f_np1 / 2 + f_n / 2) * v * dx

if problem is ProblemType.DIRICHLET:
    # apply Dirichlet boundary condition on coupling interface
    bcs.append(precice.create_coupling_dirichlet_boundary_condition(V))
if problem is ProblemType.NEUMANN:
    # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
    
    if args.method == "tr":
        F += 0.5 * precice.create_coupling_neumann_boundary_condition(v) + 0.5 * precice.create_coupling_neumann_boundary_condition(v)
    elif args.method == "ie":
        F += precice.create_coupling_neumann_boundary_condition(v)

a, L = lhs(F), rhs(F)

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
    temperature_out = File("out/%s.pvd" % precice._solver_name)
    ref_out = File("out/ref%s.pvd" % precice._solver_name)
    error_out = File("out/error%s.pvd" % precice._solver_name)
    flux_out = File("out/flux%s.pvd" % precice._solver_name)
    temperature_out << u_n
    ref_out << u_ref
    error_out << error_pointwise
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

    dt.assign(np.min([fenics_dt, precice_dt]))  # todo we could also consider deciding on time stepping size inside the adapter

    if precice_timestep_complete:
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

print(t)
u_D.t = t
u_ref = interpolate(u_D, V)
error_last, _ = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
print(error_last)
print(dt(0))

precice.finalize()

with open("errors_{participant}".format(participant=problem.name), 'a') as outfile:
    outfile.write("{}, {}\n".format(dt(0), error_last))
