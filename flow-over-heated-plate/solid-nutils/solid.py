#! /usr/bin/env python3

from nutils import cli, mesh, function, solver, export
import treelog
import numpy as np
import precice
from mpi4py import MPI


def main():

    print("Running nutils")

    # define the Nutils mesh
    grid = [np.linspace(a, b, round((b - a) / size) + 1) for (a, b, size) in [(0, 1, 0.05), (-.25, 0, 0.05)]]
    domain, geom = mesh.rectilinear(grid)

    # Nutils namespace
    ns = function.Namespace()
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
    cons = solver.optimize('lhs', sqr, droptol=1e-15)

    # preCICE setup
    participant = precice.Participant("Solid", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Solid-Mesh"
    coupling_boundary = domain.boundary['top']
    coupling_sample = coupling_boundary.sample('gauss', degree=2)  # mesh vertices at Gauss points
    vertices = coupling_sample.eval(ns.x)
    vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

    # helper functions to project heat flux to coupling boundary
    projection_matrix = coupling_boundary.integrate(ns.eval_nm('basis_n basis_m d:x'), degree=2)
    projection_cons = np.zeros(res.shape)
    projection_cons[projection_matrix.rowsupp(1e-15)] = np.nan
    def fluxdofs(v): return projection_matrix.solve(v, constrain=projection_cons)

    cons0 = cons  # to not lose the Dirichlet BC at the bottom
    lhs0 = np.zeros(res.shape)  # solution from previous timestep
    timestep = 0
    solver_dt = 0.01

    participant.initialize()

    precice_dt = participant.get_max_time_step_size()
    dt = min(precice_dt, solver_dt)

    # set u = uwall as initial condition and visualize
    sqr = domain.integral('(u - uwall)^2' @ ns, degree=2)
    lhs0 = solver.optimize('lhs', sqr)
    bezier = domain.sample('bezier', 2)
    x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs0)
    with treelog.add(treelog.DataLog()):
        export.vtk('Solid_0', bezier.tri, x, T=u)

    while participant.is_coupling_ongoing():

        # save checkpoint
        if participant.requires_writing_checkpoint():
            lhs_checkpoint = lhs0
            timestep_checkpoint = timestep

        # potentially adjust non-matching timestep sizes
        precice_dt = participant.get_max_time_step_size()
        dt = min(solver_dt, precice_dt)

        # read temperature from participant
        temperature_values = participant.read_data(mesh_name, "Temperature", vertex_ids, dt)
        temperature_function = coupling_sample.asfunction(temperature_values)

        sqr = coupling_sample.integral((ns.u - temperature_function)**2)
        cons = solver.optimize('lhs', sqr, droptol=1e-15, constrain=cons0)

        # solve nutils timestep
        lhs = solver.solve_linear('lhs', res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt))

        # write heat fluxes to participant
        flux_function = res.eval(lhs0=lhs0, lhs=lhs, dt=dt)
        flux_values = coupling_sample.eval('-flux' @ ns, fluxdofs=fluxdofs(flux_function))
        participant.write_data(mesh_name, "Heat-Flux", vertex_ids, flux_values)

        # do the coupling
        participant.advance(dt)

        # advance variables
        timestep += 1
        lhs0 = lhs

        # read checkpoint if required
        if participant.requires_reading_checkpoint():
            lhs0 = lhs_checkpoint
            timestep = timestep_checkpoint
        else:  # go to next timestep
            if timestep % 20 == 0:  # visualize
                bezier = domain.sample('bezier', 2)
                x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs)
                with treelog.add(treelog.DataLog()):
                    export.vtk('Solid_' + str(timestep), bezier.tri, x, T=u)

    participant.finalize()


if __name__ == '__main__':
    cli.run(main)
