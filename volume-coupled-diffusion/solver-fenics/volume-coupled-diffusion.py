from fenics import Function, FunctionSpace, Expression, Constant, DirichletBC, TrialFunction, TestFunction, File, \
    solve, lhs, rhs, dx, UnitSquareMesh, SubDomain, inner, grad, MeshFunction, MPI, interpolate
from fenicsprecice import Adapter
import numpy as np
import argparse


class AllDomain(SubDomain):
    def inside(self, x, on_boundary):
        return True


class AllBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] == 1


parser = argparse.ArgumentParser(description="Solving a volume coupled problem")
command_group = parser.add_mutually_exclusive_group(required=True)
command_group.add_argument("-s", "--source", help="create a source", dest="source", action="store_true")
command_group.add_argument("-d", "--drain", help="create a drain", dest="drain", action="store_true")
args = parser.parse_args()

precice = Adapter(adapter_config_filename="precice-adapter-config.json")

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, "P", 1)

u = TrialFunction(V)
v = TestFunction(V)
u_n = Function(V)
if args.source:
    u_ini = Expression("1", degree=1)
    bc = DirichletBC(V, u_ini, AllBoundary())
elif args.drain:
    u_ini = Expression("0", degree=1)
    bc = DirichletBC(V, u_ini, RightBoundary())


u_n = interpolate(u_ini, V)


precice.initialize(AllDomain(), read_function_space=V, write_object=u_n)
volume_term = precice.create_coupling_expression()
f = Function(V)

dt = precice.get_max_time_step_size()
dt_inv = Constant(1 / dt)

diffusion_source = 1
diffusion_drain = 1
if args.source:
    F = dt_inv * (u - u_n) * v * dx - (f - u_ini) * v * dx + diffusion_source * inner(grad(u), grad(v)) * dx
elif args.drain:
    F = dt_inv * (u - u_n) * v * dx - (f - u) * v * dx + diffusion_drain * inner(grad(u), grad(v)) * dx

# Time-stepping
u_np1 = Function(V)
if args.source:
    u_n.rename("Source-Data", "")
    u_np1.rename("Source-Data", "")
elif args.drain:
    u_n.rename("Drain-Data", "")
    u_np1.rename("Drain-Data", "")

t = 0

mesh_rank = MeshFunction("size_t", mesh, mesh.topology().dim())
if args.source:
    mesh_rank.set_all(MPI.rank(MPI.comm_world) + 4)
else:
    mesh_rank.set_all(MPI.rank(MPI.comm_world) + 0)
mesh_rank.rename("myRank", "")

# Generating output files
solution_out = File("output/%s.pvd" % precice.get_participant_name())
ranks = File("output/ranks%s.pvd" % precice.get_participant_name())

# output solution and reference solution at t=0, n=0
n = 0
print('output u^%d and u_ref^%d' % (n, n))
solution_out << (u_n, t)
ranks << mesh_rank

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, n)

    dt = precice.get_max_time_step_size()
    read_data = precice.read_data(dt)

    # Update the coupling expression with the new read data
    precice.update_coupling_expression(volume_term, read_data)
    f.assign(interpolate(volume_term, V))

    dt_inv.assign(1 / dt)

    # Compute solution u^n+1, use bcs u^n and coupling bcs
    a, L = lhs(F), rhs(F)
    solve(a == L, u_np1, bc)

    # Write data to preCICE according to which problem is being solved
    precice.write_data(u_np1)

    precice.advance(dt)

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
        solution_out << (u_n, t)

# Hold plot
precice.finalize()
