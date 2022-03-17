#! /usr/bin/env python3

from nutils import function, mesh, cli
import treelog as log
import numpy as np
import precice
from mpi4py import MPI


def main():

    print("Running utils")

    # define the Nutils mesh
    grid = np.linspace(0, 5, 21), np.linspace(0, 2, 9)
    domain, geom = nutils.mesh.rectilinear(grid)
    domain = domain.withboundary(inflow='left', outflow='right', wall='top,bottom') \
           - domain[8:13,:3].withboundary(wall='left,top,right')

    gauss = domain.sample('gauss', degree=2)

    # Nutils namespace
    ns = nutils.function.Namespace()
    ns.x = geom
    ns.basis = domain.basis('std', degree=1)  # linear finite elements
    ns.u = 'basis_n ?lhs_n'  # solution
    ns.dudt = 'basis_n (?lhs_n - ?lhs0_n) / ?dt'  # time derivative
    ns.vbasis = gauss.basis().vector(2)
    ns.velocity_i = 'vbasis_ni ?velocity_n'
    ns.k = 100  # thermal diffusivity
    ns.xblob = 1, 1
    ns.uinit = '.5 - .5 tanh(((x_i - xblob_i) (x_i - xblob_i) - .5) / .1)'

    # define the weak form
    res = gauss.integral('(basis_n (dudt + (velocity_i u)_,i) + k basis_n,i u_,i) d:x' @ ns)

    # define Dirichlet boundary condition
    sqr = domain.boundary['inflow'].integral('u^2 d:x' @ ns, degree=2)
    cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15)

    # preCICE setup
    interface = precice.Interface("Transport", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Transport-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    vertices = gauss.eval(ns.x)
    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)

    # coupling data
    velocity_id = interface.get_data_id("Velocity", mesh_id)

    precice_dt = interface.initialize()

    timestep = 0
    dt = 0.01

    # set u = uwall as initial condition and visualize
    sqr = domain.integral('(u - uinit)^2' @ ns, degree=2)
    lhs0 = nutils.solver.optimize('lhs', sqr)

    bezier = domain.sample('bezier', 2)
    x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs0)
    with log.add(log.DataLog()):
        nutils.export.vtk('Solid_0', bezier.tri, x, T=u)

    while interface.is_coupling_ongoing():

        # read temperature from interface
        if interface.is_read_data_available():
            velocity_values = interface.read_block_vector_data(velocity_id, vertex_ids)

        # potentially adjust non-matching timestep sizes
        dt = min(dt, precice_dt)

        # solve nutils timestep
        lhs = nutils.solver.solve_linear('lhs', res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt, velocity=velocity_values))

        # do the coupling
        precice_dt = interface.advance(dt)

        # advance variables
        timestep += 1
        lhs0 = lhs

        if timestep % 20 == 0:  # visualize
            bezier = domain.sample('bezier', 2)
            x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs)
            with log.add(log.DataLog()):
                nutils.export.vtk('Transport_' + str(timestep), bezier.tri, x, T=u)

    interface.finalize()


if __name__ == '__main__':
    nutils.cli.run(main)
