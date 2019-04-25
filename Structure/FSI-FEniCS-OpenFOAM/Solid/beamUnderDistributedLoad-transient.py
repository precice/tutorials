#Import required libs
from fenics import *
from ufl import nabla_div
import numpy as np
import matplotlib.pyplot as plt
from fenicsadapter import Adapter


#define the two kinds of boundary: clamped and coupling Neumann Boundary
def clamped_boundary(x, on_boundary):
    return on_boundary and abs(x[1])<tol

def Neumann_Boundary(x, on_boundary):
    """
    determines whether a node is on the coupling boundary
    
    """
    return on_boundary and not(abs(x[1]<tol))

# Dimensionless Geometry and material properties
d=2 #number of dimensions
H = 1
W = 0.3
rho = 0.2
E=1000.0
nu= 0.3

mu    = Constant(E / (2.0*(1.0 + nu)))

lambda_ = Constant(E*nu / ((1.0 + nu)*(1.0 - 2.0*nu)))

# create Mesh
n_x_Direction=3
n_y_Direction=10
mesh = RectangleMesh(Point(0,0), Point(W,H), n_x_Direction,n_y_Direction)


#create Function Space
V = VectorFunctionSpace(mesh, 'P', 2)

#BCs
tol=1E-14




#read fenics-adapter json-config-file
adapter_config_filename = "precice-adapter-config-fsi-s.json"

# initialize the adapter
precice = Adapter(adapter_config_filename)

precice_dt = precice.initialize(coupling_subdomain=Neumann_Boundary, 
                                mesh=mesh,
                                read_field=f_N_function,
                                write_field=u_D_function,
                                u_n=u_n)




#alpha method parameters

alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
gamma   = Constant(0.5+alpha_f-alpha_m)
beta    = Constant((gamma+0.5)**2/4.)


#parameters for Time-Stepping
T = 1.0
Nsteps = 80
dt = Constant(T/Nsteps)


#Loading by an Expression: sinus loading dependend on x_0
#p = Expression(('1*sin(2*pi*t) * x[0]','0'),degree=1, t=0)
p = Expression( ('t<1 ? 0.01 : 0.01','0'),degree=1, t=0)


# Trial and Test Functions
du = TrialFunction(V)
v = TestFunction(V)

u = Function(V)

# function known from previous timestep
u_old = Function(V)
v_old = Function(V)
a_old = Function(V)


# mark the boundary with the neumann BC
boundary_subdomains = MeshFunction("size_t", 
                                   mesh, mesh.topology().dim()-1)
boundary_subdomains.set_all(0)
force_boundary = AutoSubDomain(Neumann_Boundary)
force_boundary.mark(boundary_subdomains, 3)

#Define measure for boundary conditions integral
dss=ds(subdomain_data=boundary_subdomains)


# clamp the beam at the bottom
bc = DirichletBC(V, Constant((0,0)), clamped_boundary)


#Define strain 
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

# Define Stress tensor
def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

# Define Mass form
def m(u,v):
    return rho*inner(u,v)*dx

# Elastic stiffness form
def k(u,v):
    return inner(sigma(u), sym(grad(v))) * dx

# External Work
def Wext(u_):
    return dot(u_,p)*dss(3)


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
a_new = update_a(du, u_old, v_old, a_old, ufl=True)
v_new = update_a(a_new, u_old, v_old, a_old, ufl=True)

res = m(avg(a_old,a_new,alpha_m), v) + k(avg(u_old,du, alpha_f), v) - Wext(v)

a_form= lhs(res)
L_form= rhs(res)




# Prepare for time-stepping
time = np.linspace(0,T,Nsteps+1)
u_tip = np.zeros((Nsteps+1,1))
E_ext = 0



displacement_out = File("%s.pvd")


#time loop
for n in range(Nsteps):
    
    t=time[n+1]
    print("Time:",t)
    
    #evaluate Force field at alpha_f-average
    p.t=t-float(alpha_f*dt)
    
    #solve for new displacements
    solve(a_form==L_form, u, bc)
    
    update_fields(u, u_old, v_old, a_old)
    
    displacement_out << u
    
    p.t=t
    
    u_tip[n+1] = u(0.1,1.)[0]
    
    
    


# Plot tip displacement evolution
plt.figure()
plt.plot(time, u_tip)
plt.xlabel("Time")
plt.ylabel("Tip displacement")
plt.show()


