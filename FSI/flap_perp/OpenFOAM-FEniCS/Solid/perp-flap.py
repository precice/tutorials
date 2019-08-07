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


def extract_boundary_vertices(mesh, subdomain):
    """Extracts vertices which lie on the boundary.
    :return: stack of vertices
    """
    n = 0
    vertices_x = []
    vertices_y = []

    if not issubclass(type(subdomain), SubDomain):
        raise Exception("no correct coupling interface defined!")

    for v in dolfin.vertices(mesh):
        if subdomain.inside(v.point(), True):
            n += 1
            vertices_x.append(v.x(0))
            vertices_y.append(v.x(1))

    return np.array(vertices_x), np.array(vertices_y), n




#Specify the case you want to calculate

class StructureCase(Enum):
    OPENFOAM=1
    DUMMY2D = 2
    DUMMY3D = 3
    RFERENCE = 4
    PSREF = 5
    
Case = StructureCase.OPENFOAM

#if Case is StructureCase.OPENFOAM or StructureCase.DUMMY3D:
#    dim=2.5 #2.5 means that 2D fenics is coupled via 3d preCICE
#else:

#this is done automatically by the adapter
dim=2

#define the two kinds of boundary: clamped and coupling Neumann Boundary
def clamped_boundary(x, on_boundary):
    return on_boundary and abs(x[1])<tol

def Neumann_Boundary(x, on_boundary):
    """
    determines whether a node is on the coupling boundary
    
    """
    return on_boundary and ((abs(x[1]-1)<tol) or abs(abs(x[0])-W/2)<tol)


def top(x, on_boundary):
    return on_boundary and abs(x[1]-1)<tol
TopDomain = AutoSubDomain(top)

def left(x, on_boundary):
    return on_boundary and abs(x[0]+0.05)<tol and x[1]>=0.0
LeftDomain=AutoSubDomain(left)

# Geometry and material properties
d=2 #number of dimensions
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

if Case is not StructureCase.RFERENCE and Case is not StructureCase.PSREF:
    
    if Case is StructureCase.DUMMY2D:
        adapter_config_filename = "precice-adapter-config-fsi-dummy.json"
    elif Case is StructureCase.OPENFOAM:
        adapter_config_filename = "precice-adapter-config-fsi-s.json"
        
    elif Case is StructureCase.DUMMY3D:
        adapter_config_filename = "precice-adapter-config-fsi-dummy-3d.json"
    else: 
        raise Exception("This Case is not Implemented yet")

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

else:
    dt = Constant(0.1)



#alpha method parameters

alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
gamma   = Constant(0.5+alpha_f-alpha_m)
beta    = Constant((gamma+0.5)**2/4.)


#Loading by an Expression: sinus loading dependend on x_0
#p = Expression(('1*sin(2*pi*t) * x[0]','0'),degree=1, t=0)
#p = Expression( ('t<1 ? 0.01 : 0.01','0'),degree=1, t=0) #now precice coupling instead





force_boundary = AutoSubDomain(Neumann_Boundary)


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

if Case is StructureCase.OPENFOAM:
    Forces_x, Forces_y = precice.create_force_boundary_condition(V)

elif Case is StructureCase.DUMMY2D or Case is StructureCase.DUMMY3D:
    res -= 1/h * precice.create_coupling_neumann_boundary_condition(v)# removed the marker , 3) #3 is the marker for the domain
    #res -= dot(v, Expression( ('1','0'),degree=2)) * ds

elif Case is StructureCase.RFERENCE:
    # TODO: Maybe add boundary_markers
    

    bm= MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    bm.set_all(0)
    TopDomain.mark(bm,2)
    dss=ds(subdomain_data = bm)
    p = Constant( (F/W,0))
    res -= dot(v,p)*dss(2)

F=5


a_form= lhs(res)
L_form= rhs(res)




# Prepare for time-stepping

#parameters for Time-Stepping
#T = 1.0

t=0.0
n=0
time = []
u_tip = []
time.append(0.0)
u_tip.append(0.0)
E_ext = 0


if Case is StructureCase.DUMMY2D:
    displacement_out = File("Solid/Dummy2DOut/u_2dd.pvd")
    
elif Case is StructureCase.DUMMY3D:
    displacement_out = File("Solid/Dummy3DOut/u_3dd.pvd")
elif Case is StructureCase.OPENFOAM:
    displacement_out = File("Solid/FSI-S/u_fsi.pvd")
    
elif Case is StructureCase.RFERENCE:
    displacement_out = File("Reference/u_ref.pvd")
elif Case is StructureCase.PSREF:
    displacement_out = File("PSREF/u_psref.pvd")
    

u_n.rename("Displacement", "")
u_np1.rename("Displacement", "")
displacement_out << u_n


#time loop for coupling

if Case is StructureCase.OPENFOAM:
    
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
    
elif Case is StructureCase.DUMMY2D or Case is StructureCase.DUMMY3D:

    while precice.is_coupling_ongoing():
        
         
        #solve for new displacements
        

        solve(a_form == L_form, u_np1, bc)
        # call precice.advance
        t, n, precice_timestep_complete, precice_dt = precice.advance(u_np1, u_np1, u_n, t, float(dt), n)
        
        
        
        if precice_timestep_complete:
            
            update_fields(u_np1, saved_u_old, v_n, a_n)
            
            if n % 10==0:
                displacement_out << (u_n,t)
        
            u_tip.append(u_n(0.,1.)[0])
            time.append(t)
    
    # stop coupling    
    precice.finalize()
    
#time loop for reference calculation
elif Case is StructureCase.PSREF:
    while t<=5.5:
        
        x_coords, y_coords, n = extract_boundary_vertices(mesh, TopDomain)
        plt.scatter(x_coords, y_coords)
        A, b = assemble_system(a_form, L_form, bc)
        b_forces = b.copy()
        ps=list()
        f=2*F/n
#        for i in range(n):
#            ps.append(PointSource(V.sub(0), Point(x_coords[i], y_coords[i]), f))
#        
#        ps.append(PointSource(V.sub(0), Point(0.05,0.5),2.5))
       # ps.append(PointSource(V.sub(0), Point(0.05,1),2.5))
       # ps.append(PointSource(V.sub(0), Point(-0.05,0.5),2.5))
        x_coords, y_coords, n = extract_boundary_vertices(mesh, LeftDomain)
        
        f=2*F/n

        for i in range(n):
            ps.append(PointSource(V.sub(0), Point(x_coords[i], y_coords[i]), f))
      
        plt.scatter(x_coords, y_coords)
        for p in ps:
            p.apply(b_forces)
            
        solve(A, u_np1.vector(),b_forces)
        
        t=t+float(dt)
        update_fields(u_np1,u_n, v_n,a_n)
        displacement_out << (u_n, t)
        u_tip.append(u_n(0.,1.)[0])
        time.append(t)

elif Case is StructureCase.RFERENCE:
    while t<=5.5:
        
        solve(a_form == L_form, u_np1, bc)
            
            
        t=t+float(dt)
        update_fields(u_np1,u_n, v_n,a_n)
        displacement_out << (u_n, t)
        u_tip.append(u_n(0.,1.)[0])
        time.append(t)    

# Plot tip displacement evolution
        
displacement_out << u_n 
plt.figure()
plt.plot(time, u_tip)
plt.xlabel("Time")
plt.ylabel("Tip displacement")
plt.show()


