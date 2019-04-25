from fenics import *
from fenics import near, inner
from ufl import nabla_div, nabla_grad
import matplotlib.pyplot as plt
from fenicsadapter import *
from enum import Enum
import argparse

class ProblemType(Enum):
    """
    Enum defines problem type. Details see above.
    """
    DIRICHLET = 1  # Dirichlet problem
    NEUMANN = 2  # Neumann problem


class CouplingBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and near(x[0], x_coupling, tol):
            return True
        else:
            return False
        
def ForcesFromStresses(sigma, n_vec):
    return inner(sigma,n_vec)

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest='dirichlet', action='store_true')
parser.add_argument("-n", "--neumann", help="create a neumann problem", dest='neumann', action='store_true')

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

nx = 10
ny = 5
if problem is ProblemType.DIRICHLET:
    nx = nx*2
    adapter_config_filename = "precice-adapter-config-D.json"

elif problem is ProblemType.NEUMANN:
    adapter_config_filename = "precice-adapter-config-N.json"

#Scaled varaibles
L = 0.5 # length of the beam
W = 0.2 # width of square cross section of the beam
mu = 1 # Lame`s elasticity parameters
rho = 1 # density of the beam
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta # Lame`s elasticty parameters
g = gamma

dt = .1  # time step size
x_coupling = L  # x coordinate of coupling interface

#Mesh and function space
if problem is ProblemType.NEUMANN:
    mesh = RectangleMesh(Point(0, 0), Point(L, W), nx, ny)
elif problem is ProblemType.DIRICHLET:
    mesh = RectangleMesh(Point(L, 0), Point(L+0.5, W), nx, ny)

V = VectorFunctionSpace(mesh, 'P', 1)

#BC
tol = 1E-14
def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < tol

def coupled_boundary(x,on_boundary):
    """
    determines whether a node is on the counpling boundary
    
    """
    return on_boundary and near(x[0], x_coupling, tol)


# Define boundary condition. TODO: u_D and f_N currently dummy, just for preCICE comm. testing
u_D = Expression(("0", "0"), degree=1)
u_D_function = interpolate(u_D, V)

f_N = Expression(("x[0]", "x[1]"), degree=1)
f_N_function = interpolate(f_N, V)

boundary_markers=MeshFunction("size_t", mesh, mesh.topology().dim()-1)
coupling_boundary = CouplingBoundary()

coupling_boundary.mark(boundary_markers, 5)

#clamp the NEUMANN problem at the left end
bc_c = DirichletBC(V, Constant((0, 0)), clamped_boundary)
bcs = [bc_c]

u_n = interpolate(u_D, V) # needed for new fenicsadapter.initialize

precice = Adapter(adapter_config_filename)
if problem is ProblemType.DIRICHLET:
    precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=u_D_function,
                       write_field=f_N_function, u_n=u_n)
elif problem is ProblemType.NEUMANN:
    precice.initialize(coupling_subdomain=coupling_boundary, mesh=mesh, read_field=f_N_function,
                       write_field=u_D_function, u_n=u_n)


#Strain tensor
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    #return sym(nabla_grad(u))

def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

#Varational problem
u = TrialFunction(V)
d = u.geometric_dimension() # dimension of the space
v = TestFunction(V)
f = Constant((0, -rho*g))
T = Constant((0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds

if problem is ProblemType.DIRICHLET:
    # apply Dirichlet boundary condition on coupling interface
    bcs.append(precice.create_coupling_dirichlet_boundary_condition(V))
if problem is ProblemType.NEUMANN:
    # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
    L += precice.create_coupling_neumann_boundary_condition(v, boundary_marker=5)

u = Function(V)
# TODO: This loop is incomplete. How to handle timesteping in "static" case with preCICE? See TODO.md

#while precice.is_coupling_ongoing():
    #Solve
solve(a == L, u, bcs)
#precice.advance(u,u)

#Plot solution
plt.figure(1)
plot(u, title='Displacement')

# Plot stress
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)  # deviatoric stress
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)
plt.figure(2)
plot(von_Mises, title='Stress intensity')

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
plt.figure(3)
plot(u_magnitude, title='Displacement magnitude')
print('Displacement magnitude u min/max:',
      u_magnitude.vector().min(),
      u_magnitude.vector().max())

#Hold plot
plt.show()

#Save to file
File('output/displacement.pvd') << u
File('output/von_mises.pvd') << von_Mises
File('output/magnitude.pvd') << u_magnitude
