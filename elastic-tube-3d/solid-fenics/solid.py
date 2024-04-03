# Import required libs
from fenics import Constant, Function, AutoSubDomain, VectorFunctionSpace, interpolate, \
    TrialFunction, TestFunction, Point, Expression, DirichletBC, \
    Identity, inner, dx, ds, sym, grad, div, lhs, rhs, dot, File, solve, assemble_system
from mshr import Cylinder, generate_mesh
import numpy as np
from fenicsprecice import Adapter
import math


# define the two kinds of boundary: clamped and coupling Neumann Boundary
def clamped_boundary(x, on_boundary):
    """
    Filter nodes at both ends of tube as they are fixed
    """
    tol = 1E-14
    return on_boundary and (((abs(x[2]) - 0.0) < tol) or ((L - abs(x[2])) < tol))


def neumann_boundary(x, on_boundary):
    """
    Filter nodes which lie on the inner surface of the tube and excluding end nodes
    """
    tol = 1E-14
    return on_boundary and ((math.sqrt(x[0]**2 + x[1]**2) - R) < tol) and ((L - x[2]) > tol) and ((x[2] - 0.0) > tol)


# Geometry and material properties
dim = 3  # number of dimensions
R = 0.005
L = 0.05
rho = 1200
E = 300000
nu = 0.3

mu = Constant(E / (2.0 * (1.0 + nu)))

lambda_ = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))

# create Mesh
outer_tube = Cylinder(Point(0, 0, L), Point(0, 0, 0), R + 0.001, R + 0.001)
inner_tube = Cylinder(Point(0, 0, L), Point(0, 0, 0), R, R)
mesh = generate_mesh(outer_tube - inner_tube, 20)

# create Function Space
V = VectorFunctionSpace(mesh, 'P', 2)

# Trial and Test Functions
du = TrialFunction(V)
v = TestFunction(V)

u_np1 = Function(V)
saved_u_old = Function(V)
u_delta = Function(V)

# function known from previous timestep
u_n = Function(V)
v_n = Function(V)
a_n = Function(V)

coupling_boundary = AutoSubDomain(neumann_boundary)
fixed_boundary = AutoSubDomain(clamped_boundary)

precice = Adapter(adapter_config_filename="precice-adapter-config-fsi-s.json")

# Initialize the coupling interface
precice.initialize(coupling_boundary, read_function_space=V, write_object=V, fixed_boundary=fixed_boundary)
precice_dt = precice.get_max_time_step_size()

fenics_dt = precice_dt  # if fenics_dt == precice_dt, no subcycling is applied
dt = Constant(np.min([precice_dt, fenics_dt]))

# clamp the tube on both sides
bc = DirichletBC(V, Constant((0, 0, 0)), fixed_boundary)

# alpha method parameters
alpha_m = Constant(0)
alpha_f = Constant(0)
gamma = Constant(0.5 + alpha_f - alpha_m)
beta = Constant((gamma + 0.5) ** 2 / 4.)


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


# Functions for updating system state

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
        precice.store_checkpoint(u_n, t, n)

    precice_dt = precice.get_max_time_step_size()
    dt = Constant(np.min([precice_dt, fenics_dt]))

    # read data from preCICE and get a new coupling expression
    read_data = precice.read_data(dt)

    # Update the point sources on the coupling boundary with the new read data
    forces_x, forces_y, forces_z = precice.get_point_sources(read_data)

    A, b = assemble_system(a_form, L_form, bc)

    b_forces = b.copy()  # b is the same for every iteration, only forces change

    for ps in forces_x:
        ps.apply(b_forces)
    for ps in forces_y:
        ps.apply(b_forces)
    for ps in forces_z:
        ps.apply(b_forces)

    assert (b is not b_forces)
    solve(A, u_np1.vector(), b_forces)

    # Write relative displacements to preCICE
    u_delta.vector()[:] = u_np1.vector()[:] - u_n.vector()[:]
    precice.write_data(u_delta)

    # Call to advance coupling, also returns the optimum time step value
    precice.advance(dt(0))

    # Either revert to old step if timestep has not converged or move to next timestep
    if precice.requires_reading_checkpoint():  # roll back to checkpoint
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:
        u_n.assign(u_np1)
        t += float(dt)
        n += 1

    if precice.is_time_window_complete():
        update_fields(u_np1, saved_u_old, v_n, a_n)
        if n % 10 == 0:
            displacement_out << (u_n, t)

# Plot tip displacement evolution
displacement_out << (u_n, t)

precice.finalize()
