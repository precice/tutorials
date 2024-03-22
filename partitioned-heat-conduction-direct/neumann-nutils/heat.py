#! /usr/bin/env python3

from nutils import cli, mesh, function, solver, export
import functools
import treelog
import numpy as np
import precice


def main(n=10, degree=1, timestep=.1, alpha=3., beta=1.2):

    x_grid = np.linspace(1, 2, n)
    y_grid = np.linspace(0, 1, n)

    # define the Nutils mesh
    domain, geom = mesh.rectilinear([x_grid, y_grid])
    coupling_boundary = domain.boundary['left']
    read_sample = coupling_boundary.sample('gauss', degree=degree * 2)

    # Nutils namespace
    ns = function.Namespace()
    ns.x = geom
    ns.basis = domain.basis('std', degree=degree)
    ns.alpha = alpha  # parameter of problem
    ns.beta = beta  # parameter of problem
    ns.u = 'basis_n ?lhs_n'  # solution
    ns.dudt = 'basis_n (?lhs_n - ?lhs0_n) / ?dt'  # time derivative
    ns.flux = 'basis_n ?fluxdofs_n'  # heat flux
    ns.f = 'beta - 2 - 2 alpha'  # rhs
    ns.uexact = '1 + x_0 x_0 + alpha x_1 x_1 + beta ?t'  # analytical solution
    ns.readbasis = read_sample.basis()
    ns.readfunc = 'readbasis_n ?readdata_n'

    # define the weak form
    res = domain.integral(
        '(basis_n dudt - basis_n f + basis_n,i u_,i) d:x' @ ns, degree=degree * 2)

    # set boundary conditions at non-coupling boundaries
    # top and bottom boundary are non-coupling for both sides
    sqr = domain.boundary['top,bottom,right'].integral(
        '(u - uexact)^2 d:x' @ ns, degree=degree * 2)

    res += read_sample.integral('basis_n readfunc d:x' @ ns)

    # preCICE setup
    participant = precice.Participant("Neumann", "../precice-config.xml", 0, 1)

    mesh_name_read = "Neumann-Mesh"
    mesh_name_write = "Dirichlet-Mesh"

    vertex_ids_read = participant.set_mesh_vertices(
        mesh_name_read, read_sample.eval(ns.x))
    participant.set_mesh_access_region(mesh_name_write, [.9, 1.1, -.1, 1.1])

    participant.initialize()
    precice_dt = participant.get_max_time_step_size()
    solver_dt = timestep
    dt = min(precice_dt, solver_dt)

    vertex_ids_write, coords = participant.get_mesh_vertex_ids_and_coordinates(
        mesh_name_write)
    write_sample = domain.locate(ns.x, coords, eps=1e-10, tol=1e-10)

    precice_write = functools.partial(
        participant.write_data, mesh_name_write, "Temperature", vertex_ids_write)
    precice_read = functools.partial(
        participant.read_data, mesh_name_read, "Heat-Flux", vertex_ids_read)

    t = 0.
    istep = 0

    # initial condition
    sqr0 = domain.integral('(u - uexact)^2' @ ns, degree=degree * 2)
    lhs = solver.optimize('lhs', sqr0, arguments=dict(t=t))
    bezier = domain.sample('bezier', degree * 2)

    while participant.is_coupling_ongoing():

        # save checkpoint
        if participant.requires_writing_checkpoint():
            checkpoint = lhs, t, istep

        # prepare next timestep
        precice_dt = participant.get_max_time_step_size()
        dt = min(precice_dt, solver_dt)
        lhs0 = lhs
        istep += 1
        t += dt

        # read data from participant
        read_data = precice_read(dt)

        # update (time-dependent) boundary condition
        cons = solver.optimize('lhs', sqr, droptol=1e-15,
                               arguments=dict(t=t, readdata=read_data))

        # solve nutils timestep
        lhs = solver.solve_linear('lhs', res, constrain=cons, arguments=dict(
            lhs0=lhs0, dt=dt, t=t, readdata=read_data))

        # write data to participant
        write_data = write_sample.eval('u' @ ns, lhs=lhs)
        precice_write(write_data)

        # do the coupling
        participant.advance(dt)

        # read checkpoint if required
        if participant.requires_reading_checkpoint():
            lhs, t, istep = checkpoint
        else:
            # generate output
            x, u, uexact = bezier.eval(
                ['x_i', 'u', 'uexact'] @ ns, lhs=lhs, t=t)
            with treelog.add(treelog.DataLog()):
                export.vtk("Neumann" + "-" + str(istep), bezier.tri,
                           x, Temperature=u, reference=uexact)

    participant.finalize()


if __name__ == '__main__':
    cli.run(main)
