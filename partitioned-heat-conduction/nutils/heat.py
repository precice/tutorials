#! /usr/bin/env python3

from nutils import cli, mesh, function, solver, export
import treelog
import numpy as np
import precice


def main(side='Dirichlet'):
    print("Running nutils")

    # domain size
    y_bottom, y_top = 0, 1
    x_left, x_right = 0, 2
    x_coupling = 1  # x coordinate of coupling interface
    degree = 1  # linear finite elements
    n = 10  # number of mesh vertices per dimension

    if side == 'Dirichlet':
        x_grid = np.linspace(x_left, x_coupling, n)
    elif side == 'Neumann':
        x_grid = np.linspace(x_coupling, x_right, n)
    else:
        raise Exception('invalid side {!r}'.format(side))

    y_grid = np.linspace(y_bottom, y_top, n)

    # define the Nutils mesh
    domain, geom = mesh.rectilinear([x_grid, y_grid])
    coupling_boundary = domain.boundary['right' if side == 'Dirichlet' else 'left']
    coupling_sample = coupling_boundary.sample('gauss', degree=degree * 2)

    # Nutils namespace
    ns = function.Namespace()
    ns.x = geom
    ns.basis = domain.basis('std', degree=degree)
    ns.alpha = 3  # parameter of problem
    ns.beta = 1.3  # parameter of problem
    ns.u = 'basis_n ?lhs_n'  # solution
    ns.dudt = 'basis_n (?lhs_n - ?lhs0_n) / ?dt'  # time derivative
    ns.flux = 'basis_n ?fluxdofs_n'  # heat flux
    ns.f = 'beta - 2 - 2 alpha'  # rhs
    ns.uexact = '1 + x_0 x_0 + alpha x_1 x_1 + beta ?t'  # analytical solution
    ns.readbasis = coupling_sample.basis()
    ns.readfunc = 'readbasis_n ?readdata_n'

    # define the weak form
    res = domain.integral('(basis_n dudt - basis_n f + basis_n,i u_,i) d:x' @ ns, degree=degree * 2)

    # set boundary conditions at non-coupling boundaries
    # top and bottom boundary are non-coupling for both sides
    sqr = domain.boundary['top,bottom,left' if side == 'Dirichlet'
                     else 'top,bottom,right'].integral('(u - uexact)^2 d:x' @ ns, degree=degree * 2)

    if side == 'Dirichlet':
        sqr += coupling_sample.integral('(u - readfunc)^2 d:x' @ ns)
    else:
        res += coupling_sample.integral('basis_n readfunc d:x' @ ns)

    # preCICE setup
    interface = precice.Interface(side, "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = side + "-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    vertices = coupling_sample.eval(ns.x)
    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)

    # coupling data
    write_data = "Temperature" if side == "Neumann" else "Heat-Flux"
    read_data = "Heat-Flux" if side == "Neumann" else "Temperature"
    write_data_id = interface.get_data_id(write_data, mesh_id)
    read_data_id = interface.get_data_id(read_data, mesh_id)

    # helper functions to project heat flux to coupling boundary
    projection_matrix = coupling_boundary.integrate(ns.eval_nm('basis_n basis_m d:x'), degree=degree * 2)
    projection_cons = np.zeros(res.shape)
    projection_cons[projection_matrix.rowsupp(1e-15)] = np.nan

    def fluxdofs(v):
        return projection_matrix.solve(v, constrain=projection_cons)

    precice_dt = interface.initialize()

    # write initial data
    if interface.is_action_required(precice.action_write_initial_data()):
        write_data = np.zeros(len(vertex_ids))
        interface.write_block_scalar_data(write_data_id, vertex_ids, write_data)
        interface.mark_action_fulfilled(precice.action_write_initial_data())

    interface.initialize_data()

    t = 0

    # initial condition
    sqr0 = domain.integral('(u - uexact)^2' @ ns, degree=degree * 2)
    lhs0 = solver.optimize('lhs', sqr0, arguments=dict(t=t))
    bezier = domain.sample('bezier', degree * 2)
    x, u, uexact = bezier.eval(['x_i', 'u', 'uexact'] @ ns, lhs=lhs0, t=t)
    with treelog.add(treelog.DataLog()):
        export.vtk(side + '-0', bezier.tri, x, Temperature=u, reference=uexact)

    t += precice_dt
    timestep = 0
    dt = 0.1

    while interface.is_coupling_ongoing():

        # read data from interface
        if interface.is_read_data_available():
            readdata = interface.read_block_scalar_data(read_data_id, vertex_ids)

        # save checkpoint
        if interface.is_action_required(precice.action_write_iteration_checkpoint()):
            lhs_checkpoint = lhs0
            t_checkpoint = t
            timestep_checkpoint = timestep
            interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

        # potentially adjust non-matching timestep sizes
        dt = min(dt, precice_dt)

        # update (time-dependent) boundary condition
        cons = solver.optimize('lhs', sqr, droptol=1e-15, arguments=dict(t=t, readdata=readdata))

        # solve nutils timestep
        lhs = solver.solve_linear('lhs', res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt, t=t, readdata=readdata))

        # write data to interface
        if interface.is_write_data_required(dt):
            if side == 'Dirichlet':
                flux_function = res.eval(lhs0=lhs0, lhs=lhs, dt=dt, t=t)
                write_data = coupling_sample.eval('flux' @ ns, fluxdofs=fluxdofs(flux_function))
            else:
                write_data = coupling_sample.eval('u' @ ns, lhs=lhs)

            interface.write_block_scalar_data(write_data_id, vertex_ids, write_data)

        # do the coupling
        precice_dt = interface.advance(dt)

        # advance variables
        t += dt
        timestep += 1
        lhs0 = lhs

        # read checkpoint if required
        if interface.is_action_required(precice.action_read_iteration_checkpoint()):
            lhs0 = lhs_checkpoint
            t = t_checkpoint
            timestep = timestep_checkpoint
            interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())
        else:  # go to next timestep
            bezier = domain.sample('bezier', degree * 2)
            x, u, uexact = bezier.eval(['x_i', 'u', 'uexact'] @ ns, lhs=lhs, t=t)

            with treelog.add(treelog.DataLog()):
                export.vtk(side + "-" + str(timestep), bezier.tri, x, Temperature=u, reference=uexact)

    interface.finalize()


if __name__ == '__main__':
    cli.run(main)
