#Import required libs
from fenics import *
from fenics import Constant, Function, AutoSubDomain, RectangleMesh, VectorFunctionSpace, interpolate, \
TrialFunction, TestFunction, Point, Expression, DirichletBC, nabla_grad,\
Identity, inner,dx, ds, sym, grad, lhs, rhs, dot, File, solve, PointSource
import dolfin

from ufl import nabla_div
import numpy as np
import matplotlib.pyplot as plt
from fenicsadapter import Adapter
from enum import Enum


#define the two kinds of boundary: clamped and coupling Neumann Boundary
def clamped_boundary(x, on_boundary):
    return on_boundary and abs(x[1])<tol

def Neumann_Boundary(x, on_boundary):
    """
    determines whether a node is on the coupling boundary
    
    """
    return on_boundary and ((abs(x[1]-1)<tol) or abs(abs(x[0])-W/2)<tol)



# Geometry and material properties
dim=2 #number of dimensions
H = 1
W = 0.1
rho = 3000
E=400000.0
nu= 0.3

mu    = Constant(E / (2.0*(1.0 + nu)))

lambda_ = Constant(E*nu / ((1.0 + nu)*(1.0 - 2.0*nu)))

# create Mesh
n_x_Direction=5
n_y_Direction=50
mesh = RectangleMesh(Point(-W/2,0), Point(W/2,H), n_x_Direction,n_y_Direction)

h=Constant(H/n_y_Direction)


#create Function Space
V = VectorFunctionSpace(mesh, 'P', 2)

#BCs
tol=1E-14

# Trial and Test Functions
du = TrialFunction(V)
v = TestFunction(V)

u_np1 = Function(V)
saved_u_old = Function(V)

# function known from previous timestep
u_n = Function(V)
v_n = Function(V)
a_n = Function(V)


f_N_function = interpolate(Expression(("1","0"), degree=1), V)
u_function = interpolate(Expression(("0","0"), degree=1), V)

coupling_boundary = AutoSubDomain(Neumann_Boundary)

# create subdomain that resembles the 

## get the adapter ready

#read fenics-adapter json-config-file)

adapter_config_filename = "precice-adapter-config-fsi-s.json"
           
precice = Adapter(adapter_config_filename)

clamped_boundary_domain=AutoSubDomain(clamped_boundary)
precice_dt = precice.initialize(coupling_subdomain=coupling_boundary, 
                            mesh=mesh,
                            read_field=f_N_function,
                            write_field=u_function,
                            u_n=u_n,
                            dimension=dim,
                            dirichlet_boundary=clamped_boundary_domain)

dt = Constant(precice_dt)

force_boundary = AutoSubDomain(Neumann_Boundary)

# clamp the beam at the bottom
bc = DirichletBC(V, Constant((0,0)), clamped_boundary)

#alpha method parameters

alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
gamma   = Constant(0.5+alpha_f-alpha_m)
beta    = Constant((gamma+0.5)**2/4.)


#Define strain 
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

# Define Stress tensor
def sigma(u):
    return lambda_*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)

# Define Mass form
def m(u,v):
    return rho*inner(u,v)*dx

# Elastic stiffness form
def k(u,v):
    return inner(sigma(u), sym(grad(v))) * dx

# External Work
def Wext(u_):
    return dot(u_,p)*ds


# Update functions

# Update accelleration
def update_a(u, u_old,v_old, a_old, ufl=True):
    if ufl:
        dt_=dt
        beta_=beta
    else:
        dt_=float(dt)
        beta_=float(beta)
    
    return ((u - u_old - dt_*v_old)/beta/dt_**2 
            - (1-2*beta_)/2/beta_*a_old)

# Update velocity
def update_v(a,u_old,v_old,a_old, ufl = True):
    if ufl:
        dt_=dt
        gamma_=gamma
    else:
        dt_=float(dt)
        gamma_=float(gamma)
    
    return v_old + dt_*((1-gamma_)*a_old + gamma_*a)

def update_fields(u, u_old, v_old, a_old):
    '''Update all fields at the end of a timestep.'''
    
    u_vec, u0_vec = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()
    
    #call update functions
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    
    #assign u->u_old
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()


def avg(x_old, x_new, alpha):
    return alpha*x_old + (1-alpha)*x_new

# residual
a_np1 = update_a(du, u_n, v_n, a_n, ufl=True)
v_np1 = update_v(a_np1, u_n, v_n, a_n, ufl=True)

res = m(avg(a_n,a_np1,alpha_m), v) + k(avg(u_n,du, alpha_f), v)  #TODO: Wext(v) needs to be replaced by coupling

Forces_x, Forces_y = precice.create_force_boundary_condition(V)

a_form= lhs(res)
L_form= rhs(res)

# Prepare for time-stepping

#parameters for Time-Stepping
t=0.0
n=0
time = []
u_tip = []
time.append(0.0)
u_tip.append(0.0)
E_ext = 0


displacement_out = File("Solid/FSI-S/u_fsi.pvd")
    
u_n.rename("Displacement", "")
u_np1.rename("Displacement", "")
displacement_out << u_n


#time loop for coupling

  
while precice.is_coupling_ongoing():
    A, b = assemble_system(a_form, L_form, bc)

    b_forces = b.copy() # b is the same for every iteration, only forces change
    
    for ps in Forces_x:
        ps.apply(b_forces)
    for ps in Forces_y:
        ps.apply(b_forces)
        
    assert(b is not b_forces)
    solve(A, u_np1.vector(), b_forces)
    
    t, n, precice_timestep_complete, precice_dt, Forces_x, Forces_y = precice.advance(u_np1, u_np1, u_n, t, float(dt), n)
    
    
    if precice_timestep_complete:
        
        update_fields(u_np1, saved_u_old, v_n, a_n)
        
        if n % 10==0:
            displacement_out << (u_n,t)
    
        u_tip.append(u_n(0.,1.)[0])
        time.append(t)
    
# Plot tip displacement evolution
        
displacement_out << u_n 
plt.figure()
plt.plot(time, u_tip)
plt.xlabel("Time")
plt.ylabel("Tip displacement")
plt.show()


