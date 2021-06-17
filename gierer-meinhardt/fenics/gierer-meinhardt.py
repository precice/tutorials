"""
This basic example implements the Activator in the Gierer-Meinhardt model following http://www.scholarpedia.org/article/Gierer-Meinhardt_model

The code is implemented following https://odr.chalmers.se/bitstream/20.500.12380/300668/1/Sebastian_Persson_Master_thesis.pdf and https://github.com/sebapersson/Master_thesis, https://github.com/sebapersson/Master_thesis/blob/1643193a8f786674ce0573518162c0497306935d/Code/Python/RD_models/Simulate_RD_models.py#L289-L300
"""

from fenics import Function, FunctionSpace, Expression, Constant, DirichletBC, TrialFunction, TestFunction, \
    File, solve, lhs, rhs, dx, UnitSquareMesh, SubDomain, inner, grad, MeshFunction, MPI, interpolate, assemble, derivative, FiniteElement, MixedElement, dot
from fenicsprecice import Adapter
import numpy as np
import argparse


class AllDomain(SubDomain):
    def inside(self, x, on_boundary):
        return True


parser = argparse.ArgumentParser(description="Solving heat equation for simple or complex interface case")
command_group = parser.add_mutually_exclusive_group(required=True)
command_group.add_argument("-a", "--activator", help="create a activator", dest="activator",
                           action="store_true")
command_group.add_argument("-i", "--inhibitor", help="create a inhibitor", dest="inhibitor", action="store_true")
args = parser.parse_args()

if args.activator:
    precice = Adapter(adapter_config_filename="precice-adapter-config-activator.json")
elif args.inhibitor:
    precice = Adapter(adapter_config_filename="precice-adapter-config-inhibitor.json")

# from https://github.com/sebapersson/Master_thesis/blob/1643193a8f786674ce0573518162c0497306935d/Code/Python/RD_models/Run_RD_simulations.py#L21
a = Constant(0.5)
b = Constant(2.0)
gamma = Constant(20)
d = Constant(50)

mesh = UnitSquareMesh(10, 10)

# Standard Lagrange elements
#P1 = FiniteElement('P', mesh.ufl_cell(), 1)
#element = MixedElement([P1, P1])
#W = FunctionSpace(mesh, element)
#u_mixed_n = Function(W)
#u_n, f = u_mixed_n.split()

V = FunctionSpace(mesh, "P", 1)

#u = Function(V)  # see https://fenicsproject.discourse.group/t/adding-expressions-with-non-matching-form-arguments-vs-v-1/1867
u = TrialFunction(V)

v = TestFunction(V)
u_n = Function(V)
#u_ini = Expression("x[0]*x[1]", degree=1)
#u_n = interpolate(u_ini, V)

dt = precice.initialize(AllDomain(), read_function_space=V, write_object=u_n)
volume_term = Expression("x[0]*x[1]", degree=1) 
#volume_term = precice.create_coupling_expression()

print(precice._fenics_vertices)

dt_inv = Constant(1/dt)

if args.activator:
    F = dt_inv*(u - u_n)*v*dx - gamma*(a - u + u*volume_term)*v*dx + inner(grad(u), grad(v))*dx
    #F = dt_inv*(u - u_n)*v*dx - gamma*(a - u + u**2*volume_term)*v*dx + inner(grad(u), grad(v))*dx
elif args.inhibitor:
    F = dt_inv*(u - u_n)*v*dx - gamma*(b - (volume_term**2)*u)*v*dx + d*inner(grad(u), grad(v))*dx

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
if args.activator:
    u_np1.rename("Substance", "")
elif args.inhibitor:
    u_np1.rename("Antagonist", "")

t = 0

mesh_rank = MeshFunction("size_t", mesh, mesh.topology().dim())
if args.activator:
    mesh_rank.set_all(MPI.rank(MPI.comm_world) + 4)
else:
    mesh_rank.set_all(MPI.rank(MPI.comm_world) + 0)
mesh_rank.rename("myRank", "")

# Generating output files
solution_out = File("out/%s.pvd" % precice.get_participant_name())
ranks = File("out/ranks%s.pvd" % precice.get_participant_name())

# output solution and reference solution at t=0, n=0
n = 0
print('output u^%d and u_ref^%d' % (n, n))
solution_out << u_n
ranks << mesh_rank

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.is_action_required(precice.action_write_iteration_checkpoint()):
        precice.store_checkpoint(u_n, t, n)

    read_data = precice.read_data()    

    # Update the coupling expression with the new read data
    #precice.update_coupling_expression(volume_term, read_data)

    dt_inv.assign(1/dt)

    # Compute solution u^n+1, use bcs u^n and coupling bcs
    solve(a == L, u_np1)
    #solve(F == 0, u_np1)

    """
    # From https://github.com/sebapersson/Master_thesis/blob/1643193a8f786674ce0573518162c0497306935d/Code/Python/RD_models/Simulate_RD_models.py#L388-L395
    # Iterate 10-times when trying to solve the problem
    # In the future look into explicit stepping 
    for m in range(0, 10):
        A = assemble(derivative(F, u), keep_diagonal = True) # Define the derivative
        R = assemble(F) # Define the right hand side
        A.ident_zeros() # Fix "zero"-rows to get consistent system
        solve(A, dU.vector(), R) # Solver matrix system
        u.assign(u - dU) # Update solution one step  
    """

    # Write data to preCICE according to which problem is being solved
    precice.write_data(u_np1)

    dt = precice.advance(dt)

    # roll back to checkpoint
    if precice.is_action_required(precice.action_read_iteration_checkpoint()):
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:  # update solution
        u_n.assign(u_np1)
        t += float(dt)
        n += 1

    if precice.is_time_window_complete():
        solution_out << u_n

# Hold plot
precice.finalize()
