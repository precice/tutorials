#!/usr/bin/env python
# coding: utf-8

# The basic example is taken from *Langtangen, Hans Petter, and Anders Logg. Solving PDEs in Python: The FEniCS
# Tutorial I. Springer International Publishing, 2016.*
#
# The example code has been extended with preCICE API calls to allow for an alternating Schwarz  overlapping coupling of two separate heat equations.
#
# The original source code can be found on https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft03_heat.py.
#
# Heat equation with Dirichlet conditions. (Dirichlet problem, $f_N = \mathcal{D}\left(u_C\right)$)
#
# \begin{align*}
# \frac{\partial u}{\partial t} &= \Delta u + f && \text{on the domain $\Omega = \left[0,1+\delta\right] \times \left[0,1\right]$}\\
# u &= u_C             && \text{on the coupling boundary $\Gamma_C$ at $x = 1+\delta$}\\
# u &= u_D             && \text{on the remaining boundary $\Gamma = \partial \Omega \setminus \Gamma_C$}\\
# u &= u_0             && \text{on $\Omega$ at t = 0}\\
# u &= 1 + x^2 + \alpha y^2 + \beta t \\
# f &= \beta - 2 - 2\alpha \\
# \end{align*}
#
# Similarly a second heat equation is defined on $\Omega = \left[1-\delta,2\right]$ with dirichlet boundary conditions applied at $1-\delta$.
# The alternating Schwarz method is applied in the overlap domain where $x \in \left[1-\delta, 1+\delta\right]$.

from __future__ import print_function, division
from fenics import Function, FunctionSpace, Expression, Constant, DirichletBC, TrialFunction, TestFunction, File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate
from fenics import SubDomain, Point, RectangleMesh, near, Function, Expression
from fenicsprecice import Adapter
from errorcomputation import compute_errors
import argparse
import numpy as np
from dolfin import dot
from enum import Enum
import sys


class Participant(Enum):
    """
    Enum defines which part of the domain [x_left, x_right] x [y_bottom, y_top] we compute.
    """
    LEFT = "Left"  # left part of domain in simple interface case
    RIGHT = "Right"  # right part of domain in simple interface case


parser = argparse.ArgumentParser(description="Solving heat equation with Schwarz type domain decomposition")
parser.add_argument("participantName", help="Name of the solver.", type=str, choices=[p.value for p in Participant])
args = parser.parse_args()

fenics_dt = .1  # time step size
# Error is bounded by coupling accuracy. In theory we would obtain the analytical solution.
error_tol = 10**-6
alpha = 3  # parameter alpha
beta = 1.3  # parameter beta

participant_name = args.participantName

y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
# x coordinate of coupling interface; for Schwarz Domain Decomposition coupling interface sits between The DoFs.
x_coupling = 1.0

overlap_cells = 1

nx = 9 + overlap_cells
ny = 9
hx = (x_right - x_left) / (2 * nx - overlap_cells)


if participant_name == Participant.LEFT.value:
    p0 = Point(x_left, y_bottom)
    # rightmost point is the DoF hx/2 right of the coupling interface that belongs to the Right participant
    x_schwarz_read = x_coupling + hx * (overlap_cells * 0.5)
    x_schwarz_write = x_coupling - hx * (overlap_cells * 0.5)
    p1 = Point(x_schwarz_read, y_top)
elif participant_name == Participant.RIGHT.value:
    # leftmost point is the DoF hx/2 left of the coupling interface that belongs to the Left participant
    x_schwarz_read = x_coupling - hx * (overlap_cells * 0.5)
    x_schwarz_write = x_coupling + hx * (overlap_cells * 0.5)
    p0 = Point(x_schwarz_read, y_bottom)
    p1 = Point(x_right, y_top)


class ExcludeStraightBoundary(SubDomain):
    def get_user_input_args(self, args):
        self._interface = args.interface

    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and not near(x[0], x_schwarz_read, tol) or near(x[1], y_top, tol) or near(x[1], y_bottom, tol):
            return True
        else:
            return False


class OverlapDomain(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if (x[0] <= x_coupling + hx * (overlap_cells * 0.5) + tol) and (x[0] >= x_coupling -
                                                                        hx * (overlap_cells - 0.5) - tol):  # Point lies inside of overlapping domain
            return True
        else:
            return False


class ReadBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and near(x[0], x_schwarz_read, tol):
            return True
        else:
            return False


class WriteBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if near(x[0], x_schwarz_write, tol):  # Point lies inside of the domain!
            return True
        else:
            return False


mesh = RectangleMesh(p0, p1, nx, ny, diagonal="left")
read_boundary = ReadBoundary()
write_boundary = WriteBoundary()
remaining_boundary = ExcludeStraightBoundary()


# Define function space using mesh
V = FunctionSpace(mesh, 'P', 2)

# Define boundary conditions
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0)
u_D_function = interpolate(u_D, V)

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")


precice, precice_dt, initial_data = None, 0.0, None

# Initialize the adapter according to the specific participant
precice = Adapter(adapter_config_filename="precice-adapter-config.json")

precice.initialize(OverlapDomain(), read_function_space=V, write_object=u_D_function)

precice_dt = precice.get_max_time_step_size()
dt = Constant(0)
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('beta - 2 - 2*alpha', degree=2, alpha=alpha, beta=beta, t=0)
F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f) * v * dx

bcs = [DirichletBC(V, u_D, remaining_boundary)]
coupling_expression = precice.create_coupling_expression()
bcs.append(DirichletBC(V, coupling_expression, read_boundary))

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
u_np1.rename("Temperature", "")
t = 0

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

# Generating output files
temperature_out = File("output/%s.pvd" % precice.get_participant_name())
ref_out = File("output/ref%s.pvd" % precice.get_participant_name())
error_out = File("output/error%s.pvd" % precice.get_participant_name())

# output solution and reference solution at t=0, n=0
n = 0
print('output u^%d and u_ref^%d' % (n, n))
temperature_out << u_n
ref_out << u_ref

error_total, error_pointwise = compute_errors(u_n, u_ref, V)
error_out << error_pointwise

# set t_1 = t_0 + dt, this gives u_D^1
# call dt(0) to evaluate FEniCS Constant. Todo: is there a better way?
u_D.t = t + dt(0)
f.t = t + dt(0)

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, n)

    precice_dt = precice.get_max_time_step_size()
    dt.assign(np.min([fenics_dt, precice_dt]))

    read_data = precice.read_data(dt)

    # Update the coupling expression with the new read data
    precice.update_coupling_expression(coupling_expression, read_data)

    # Compute solution u^n+1, use bcs u_D^n+1, u^n and coupling bcs
    solve(a == L, u_np1, bcs)

    precice.write_data(u_np1)

    precice_dt = precice.advance(dt(0))

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

    temperature_out << u_np1
    ref_out << u_ref
    error_out << error_pointwise

    # Update Dirichlet BC
    u_D.t = t + float(dt)
    f.t = t + float(dt)

# Hold plot
precice.finalize()
