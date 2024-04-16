"""
This source code is based on https://gitlab.enpc.fr/jeremy.bleyer/comet-fenics/-/blob/master/examples/elastodynamics/demo_elastodynamics.py.rst

See also:
    Jeremy Bleyer. (2018). Numerical Tours of Computational Mechanics with FEniCS. Zenodo. https://doi.org/10.5281/zenodo.1287832
"""

# Import required libs
from fenics import Constant, Function, AutoSubDomain, RectangleMesh, VectorFunctionSpace, interpolate, \
    TrialFunction, TestFunction, Point, Expression, DirichletBC, project, \
    Identity, inner, dx, ds, sym, grad, div, lhs, rhs, dot, File, solve, assemble_system
import numpy as np
import matplotlib.pyplot as plt
from fenicsprecice import Adapter
from enum import Enum


# define the two kinds of boundary: clamped and coupling Neumann Boundary
def clamped_boundary(x, on_boundary):
    return on_boundary and abs(x[1]) < tol


def neumann_boundary(x, on_boundary):
    """
    determines whether a node is on the coupling boundary

    """
    return on_boundary and ((abs(x[1] - 1) < tol) or abs(abs(x[0]) - W / 2) < tol)


# Geometry and material properties
dim = 2  # number of dimensions
H = 1
W = 0.1
rho = 3000
E = 4000000
nu = 0.3

mu = Constant(E / (2.0 * (1.0 + nu)))

lambda_ = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))

# create Mesh
n_x_Direction = 4
n_y_Direction = 26
mesh = RectangleMesh(Point(-W / 2, 0), Point(W / 2, H), n_x_Direction, n_y_Direction)

h = Constant(H / n_y_Direction)

# create Function Space
V = VectorFunctionSpace(mesh, 'P', 2)

# BCs
tol = 1E-14

# Trial and Test Functions
du = TrialFunction(V)
v = TestFunction(V)

u_np1 = Function(V)

# function known from previous timestep
u_n = Function(V)
v_n = Function(V)
a_n = Function(V)


f_N_function = interpolate(Expression(("1", "0"), degree=1), V)
u_function = interpolate(Expression(("0", "0"), degree=1), V)

coupling_boundary = AutoSubDomain(neumann_boundary)
fixed_boundary = AutoSubDomain(clamped_boundary)

precice = Adapter(adapter_config_filename="precice-adapter-config-fsi-s.json")

# Initialize the coupling interface
precice.initialize(coupling_boundary, read_function_space=V, write_object=V, fixed_boundary=fixed_boundary)

precice_dt = precice.get_max_time_step_size()
fenics_dt = precice_dt  # if fenics_dt == precice_dt, no subcycling is applied
# n_substeps = 5  # number of substeps per window
# fenics_dt = precice_dt / n_substeps  # if fenics_dt < precice_dt, subcycling is applied
dt = Constant(np.min([precice_dt, fenics_dt]))

# clamp the beam at the bottom
bc = DirichletBC(V, Constant((0, 0)), fixed_boundary)

# alpha method parameters
alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
# alpha_m = Constant(0)
# alpha_f = Constant(0)

"""
Check requirements for alpha_m and alpha_f from
    Chung, J., and Hulbert, G. M. (June 1, 1993). "A Time Integration Algorithm for Structural Dynamics With Improved Numerical Dissipation:
    The Generalized-α Method." ASME. J. Appl. Mech. June 1993; 60(2): 371–375. https://doi.org/10.1115/1.2900803
"""
assert (float(alpha_m) <= float(alpha_f))
assert (float(alpha_f) <= 0.5)

gamma = Constant(0.5 + alpha_f - alpha_m)
beta = Constant((gamma + 0.5)**2 / 4.)


# Define strain
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)


# Define Stress tensor
def sigma(u):
    return lambda_ * div(u) * Identity(dim) + 2 * mu * epsilon(u)


# Define Mass form
def m(u, v):
    return rho * inner(u, v) * dx


# Elastic stiffness form
def k(u, v):
    return inner(sigma(u), sym(grad(v))) * dx


# External Work
def Wext(u_):
    return dot(u_, p) * ds


# Update functions

# Update acceleration
def update_a(u, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)

    return ((u - u_old - dt_ * v_old) / beta / dt_ ** 2
            - (1 - 2 * beta_) / 2 / beta_ * a_old)


# Update velocity
def update_v(a, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)

    return v_old + dt_ * ((1 - gamma_) * a_old + gamma_ * a)


def update_fields(u, u_old, v_old, a_old):
    """Update all fields at the end of a timestep."""

    u_vec, u0_vec = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()

    # call update functions
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # assign u->u_old
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()


def avg(x_old, x_new, alpha):
    return alpha * x_old + (1 - alpha) * x_new


# residual
a_np1 = update_a(du, u_n, v_n, a_n, ufl=True)
v_np1 = update_v(a_np1, u_n, v_n, a_n, ufl=True)

res = m(avg(a_n, a_np1, alpha_m), v) + k(avg(u_n, du, alpha_f), v)

a_form = lhs(res)
L_form = rhs(res)

# parameters for Time-Stepping
t = 0.0
n = 0
E_ext = 0

displacement_out = File("output/u_fsi.pvd")

u_n.rename("Displacement", "")
u_np1.rename("Displacement", "")
displacement_out << (u_n, t)

while precice.is_coupling_ongoing():

    if precice.requires_writing_checkpoint():  # write checkpoint
        precice.store_checkpoint((u_n, v_n, a_n), t, n)

    precice_dt = precice.get_max_time_step_size()
    dt = Constant(np.min([precice_dt, fenics_dt]))

    # read data from preCICE and get a new coupling expression
    # sample force F at $F(t_{n+1-\alpha_f})$ (see generalized alpha paper)
    read_data = precice.read_data((1 - float(alpha_f)) * dt)

    # Update the point sources on the coupling boundary with the new read data
    Forces_x, Forces_y = precice.get_point_sources(read_data)

    A, b = assemble_system(a_form, L_form, bc)

    b_forces = b.copy()  # b is the same for every iteration, only forces change

    for ps in Forces_x:
        ps.apply(b_forces)
    for ps in Forces_y:
        ps.apply(b_forces)

    assert (b is not b_forces)
    solve(A, u_np1.vector(), b_forces)

    # Write new displacements to preCICE
    precice.write_data(u_np1)

    # Call to advance coupling, also returns the optimum time step value
    precice.advance(float(dt))

    # Either revert to old step if timestep has not converged or move to next timestep
    if precice.requires_reading_checkpoint():  # roll back to checkpoint
        uva_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_cp, v_cp, a_cp = uva_cp
        u_n.assign(u_cp)
        v_n.assign(v_cp)
        a_n.assign(a_cp)
        t = t_cp
        n = n_cp
    else:
        update_fields(u_np1, u_n, v_n, a_n)
        u_n.assign(u_np1)
        t += float(dt)
        n += 1

    if precice.is_time_window_complete():
        if n % 10 == 0:
            displacement_out << (u_n, t)

# Plot tip displacement evolution
displacement_out << (u_n, t)

precice.finalize()
