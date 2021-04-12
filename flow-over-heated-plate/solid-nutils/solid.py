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
    interface = precice.Interface("Solid", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Solid-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    coupling_boundary = domain.boundary['top']
    coupling_sample = coupling_boundary.sample('gauss', degree=2)  # mesh vertices at Gauss points
    vertices = coupling_sample.eval(ns.x)
    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)

    # coupling data
    flux_id = interface.get_data_id("Heat-Flux", mesh_id)
    temperature_id = interface.get_data_id("Temperature", mesh_id)

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
            temperature_values = interface.read_block_scalar_data(temperature_id, vertex_ids)
            temperature_function = coupling_sample.asfunction(temperature_values)

            sqr = coupling_sample.integral((ns.u - temperature_function)**2)
            cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15, constrain=cons0)

        # save checkpoint
        if interface.is_action_required(precice.action_write_iteration_checkpoint()):
            lhs_checkpoint = lhs0
            timestep_checkpoint = timestep
            interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

        # potentially adjust non-matching timestep sizes
        dt = min(dt, precice_dt)

        # solve nutils timestep
        lhs = nutils.solver.solve_linear('lhs', res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt))

        # write heat fluxes to interface
        if interface.is_write_data_required(dt):
            flux_function = res.eval(lhs0=lhs0, lhs=lhs, dt=dt)
            flux_values = coupling_sample.eval('-flux' @ ns, fluxdofs=fluxdofs(flux_function))
            interface.write_block_scalar_data(flux_id, vertex_ids, flux_values)

        # do the coupling
        precice_dt = interface.advance(dt)

        # advance variables
        timestep += 1
        lhs0 = lhs

        # read checkpoint if required
        if interface.is_action_required(precice.action_read_iteration_checkpoint()):
            lhs0 = lhs_checkpoint
            timestep = timestep_checkpoint
            interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())
        else:  # go to next timestep
            if timestep % 20 == 0:  # visualize
                bezier = domain.sample('bezier', 2)
                x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs)
                with treelog.add(treelog.DataLog()):
                    nutils.export.vtk('Solid_' + str(timestep), bezier.tri, x, T=u)

    interface.finalize()


if __name__ == '__main__':
    nutils.cli.run(main)
