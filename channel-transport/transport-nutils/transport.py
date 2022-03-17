#! /usr/bin/env python3

import nutils
import treelog
import numpy as np
import precice
from mpi4py import MPI


def main():

    print("Running utils")

    # define the Nutils mesh
    grid = [np.linspace(a, b, round((b - a) / size) + 1) for (a, b, size) in [(0, 1, 0.05), (-.25, 0, 0.05)]]
    domain, geom = nutils.mesh.rectilinear(grid)

    # Nutils namespace
    ns = nutils.function.Namespace()
    ns.x = geom

    ns.basis = domain.basis('std', degree=1)  # linear finite elements
    ns.u = 'basis_n ?lhs_n'  # solution
    ns.dudt = 'basis_n (?lhs_n - ?lhs0_n) / ?dt'  # time derivative
    ns.flux = 'basis_n ?fluxdofs_n'  # heat flux
    ns.k = 100  # thermal diffusivity
    ns.uwall = 310  # wall temperature

    # define the weak form
    res = domain.integral('(basis_n dudt + k basis_n,i u_,i) d:x' @ ns, degree=2)

    # define Dirichlet boundary condition
    sqr = domain.boundary['bottom'].integral('(u - uwall)^2 d:x' @ ns, degree=2)
    cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15)

    # preCICE setup
    interface = precice.Interface("Transport", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Transport-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    coupling_boundary = domain.boundary['top']
    coupling_sample = coupling_boundary.sample('gauss', degree=2)  # mesh vertices at Gauss points
    vertices = coupling_sample.eval(ns.x)
    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)

    # coupling data
    velocity_id = interface.get_data_id("Velocity", mesh_id)

    # helper functions to project heat flux to coupling boundary
    projection_matrix = coupling_boundary.integrate(ns.eval_nm('basis_n basis_m d:x'), degree=2)
    projection_cons = np.zeros(res.shape)
    projection_cons[projection_matrix.rowsupp(1e-15)] = np.nan
    def fluxdofs(v): return projection_matrix.solve(v, constrain=projection_cons)

    precice_dt = interface.initialize()

    cons0 = cons  # to not lose the Dirichlet BC at the bottom
    lhs0 = np.zeros(res.shape)  # solution from previous timestep
    timestep = 0
    dt = 0.01

    # set u = uwall as initial condition and visualize
    sqr = domain.integral('(u - uwall)^2' @ ns, degree=2)
    lhs0 = nutils.solver.optimize('lhs', sqr)
    bezier = domain.sample('bezier', 2)
    x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs0)
    with treelog.add(treelog.DataLog()):
        nutils.export.vtk('Solid_0', bezier.tri, x, T=u)

    while interface.is_coupling_ongoing():

        # read temperature from interface
        if interface.is_read_data_available():
            velocity_values = interface.read_block_vector_data(velocity_id, vertex_ids)
            velocity_function = coupling_sample.asfunction(velocity_values)

            sqr = coupling_sample.integral((ns.u - velocity_function)**2)
            cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15, constrain=cons0)

        # potentially adjust non-matching timestep sizes
        dt = min(dt, precice_dt)

        # solve nutils timestep
        lhs = nutils.solver.solve_linear('lhs', res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt))

        # do the coupling
        precice_dt = interface.advance(dt)

        # advance variables
        timestep += 1
        lhs0 = lhs

        if timestep % 20 == 0:  # visualize
            bezier = domain.sample('bezier', 2)
            x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs)
            with treelog.add(treelog.DataLog()):
                nutils.export.vtk('Transport_' + str(timestep), bezier.tri, x, T=u)

    interface.finalize()


if __name__ == '__main__':
    nutils.cli.run(main)
