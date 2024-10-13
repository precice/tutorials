#!/usr/bin/env python3

from nutils import mesh, function, solver, cli
from nutils.expression_v2 import Namespace
import numpy
import treelog
import matplotlib.pyplot as plt
import precice


def main(nelems: int, etype: str, degree: int, reynolds: float):
    '''
    1D channel flow problem.
    .. arguments::
       nelems [12]
         Number of elements along edge.
       etype [square]
         Element type (square/triangle/mixed).
       degree [2]
         Polynomial degree for velocity; the pressure space is one degree less.
       reynolds [1000]
         Reynolds number, taking the domain size as characteristic length.
    '''
    # preCICE setup
    participant_name = "Fluid1D"
    config_file_name = "../precice-config.xml"
    solver_process_index = 0
    solver_process_size = 1
    interface = precice.Interface(participant_name, config_file_name, solver_process_index, solver_process_size)
    mesh_name = "Fluid1D-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    velocity_name = "Velocity"
    velocity_id = interface.get_data_id(velocity_name, mesh_id)
    pressure_name = "Pressure"
    pressure_id = interface.get_data_id(pressure_name, mesh_id)
    positions = [0, 0, 0]
    vertex_ids = interface.set_mesh_vertex(mesh_id, positions)
    precice_dt = interface.initialize()

    # problem definition
    domain, geom = mesh.rectilinear([numpy.linspace(0, 1, nelems + 1)])
    ns = Namespace()
    ns.δ = function.eye(domain.ndims)
    ns.Σ = function.ones([domain.ndims])
    ns.Re = reynolds
    ns.x = geom
    ns.define_for('x', gradient='∇', normal='n', jacobians=('dV', 'dS'))
    ns.ubasis = domain.basis('std', degree=2).vector(domain.ndims)
    ns.pbasis = domain.basis('std', degree=1)
    ns.u = function.dotarg('u', ns.ubasis)
    ns.p = function.dotarg('p', ns.pbasis)
    ns.stress_ij = '(∇_j(u_i) + ∇_i(u_j)) / Re - p δ_ij'
    ures = domain.integral('∇_j(ubasis_ni) stress_ij dV' @ ns, degree=4)
    pres = domain.integral('pbasis_n ∇_k(u_k) dV' @ ns, degree=4)

    while interface.is_coupling_ongoing():
        if interface.is_read_data_available():  # get dirichlet pressure outlet value from 3D solver
            p_read = interface.read_scalar_data(pressure_id, vertex_ids)
            p_read = numpy.maximum(0, p_read)     # filter out unphysical negative pressure values
        else:
            p_read = 0
        usqr = domain.boundary['left'].integral('(u_0 - 1)^2 dS' @ ns, degree=2)
        diricons = solver.optimize('u', usqr, droptol=1e-15)
        ucons = diricons
        stringintegral = '(p - ' + str(numpy.rint(p_read)) + ')^2 dS'
        psqr = domain.boundary['right'].integral(stringintegral @ ns, degree=2)
        pcons = solver.optimize('p', psqr, droptol=1e-15)
        cons = dict(u=ucons, p=pcons)
        with treelog.context('stokes'):
            state0 = solver.solve_linear(('u', 'p'), (ures, pres), constrain=cons)
            x, u, p = postprocess(domain, ns, precice_dt, **state0)

        if interface.is_action_required(
                precice.action_write_iteration_checkpoint()):
            u_iter = u
            p_iter = p
            interface.mark_action_fulfilled(
                precice.action_write_iteration_checkpoint())

        write_vel = [0, 0, u[-1]]
        if interface.is_write_data_required(precice_dt):    # write new velocities to 3D solver
            interface.write_vector_data(velocity_id, vertex_ids, write_vel)
        precice_dt = interface.advance(precice_dt)

        if interface.is_action_required(
                precice.action_read_iteration_checkpoint()):
            u = u_iter
            p = p_iter
            interface.mark_action_fulfilled(
                precice.action_read_iteration_checkpoint())

    interface.finalize()
    return state0


def postprocess(domain, ns, dt, **arguments):

    bezier = domain.sample('bezier', 9)
    x, u, p = bezier.eval(['x_i', 'u_i', 'p'] @ ns, **arguments)

    ax1 = plt.subplot(211)
    ax1.plot(x, u)
    ax1.set_title("velocity u")
    ax2 = plt.subplot(212)
    ax2.plot(x, p)
    ax2.set_title("pressure p")
    plt.savefig("./results/Fluid1D_" + str(dt) + ".png")
    return x, u, p


if __name__ == '__main__':
    cli.run(main)
