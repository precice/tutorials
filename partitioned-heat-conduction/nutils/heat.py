#! /usr/bin/env python3

import nutils
import numpy
import treelog
import precice
from precice import action_write_initial_data, action_read_iteration_checkpoint, action_write_iteration_checkpoint


def main(nelems: 'number of elements along edge' = 10,
         btype: 'type of basis function (std/spline)' = 'std',
         degree: 'polynomial degree' = 1,
         side='Dirichlet',
         timestepsize=.1):

    print("Running nutils")

    y_bottom, y_top = 0, 1
    x_left, x_right = 0, 2
    x_coupling = 1  # x coordinate of coupling interface

    if side == 'Dirichlet':
        x_grid = numpy.linspace(x_left, x_coupling, nelems)
    elif side == 'Neumann':
        x_grid = numpy.linspace(x_coupling, x_right, nelems)
    else:
        raise Exception('invalid side {!r}'.format(side))

    y_grid = numpy.linspace(y_bottom, y_top, nelems)

    domain, geom = nutils.mesh.rectilinear([x_grid, y_grid])

    ns = nutils.function.Namespace()
    ns.x = geom
    ns.basis = domain.basis(btype, degree=degree)
    ns.alpha = 3
    ns.beta = 1.3
    ns.u = 'basis_n ?lhs_n'
    ns.dudt = 'basis_n (?lhs_n - ?lhs0_n) / ?dt'
    ns.flux = 'basis_n ?fluxdofs_n'
    ns.f = 'beta - 2 - 2 alpha'
    ns.uexact = '1 + x_0 x_0 + alpha x_1 x_1 + beta ?t'

    res0 = domain.integral(
        '(basis_n dudt - basis_n f + basis_n,i u_,i) d:x' @ ns, degree=degree * 2)

    # set boundary conditions at non-coupling boundaries
    # top and bottom boundary are non-coupling for both sides
    sqr0 = domain.boundary['top'].integral(
        '(u - 1 - x_0 x_0 - alpha - beta ?t)^2 d:x' @ ns, degree=degree * 2)
    sqr0 += domain.boundary['bottom'].integral(
        '(u - 1 - x_0 x_0 - beta ?t)^2 d:x' @ ns, degree=degree * 2)
    if side == 'Dirichlet':  # left boundary is non-coupling
        sqr0 += domain.boundary['left'].integral(
            '(u - 1 - alpha x_1 x_1 - beta ?t)^2 d:x' @ ns, degree=degree * 2)
    elif side == 'Neumann':  # right boundary is non-coupling
        sqr0 += domain.boundary['right'].integral(
            '(u - 1 - x_0 x_0 - alpha x_1 x_1 - beta ?t)^2 d:x' @ ns, degree=degree * 2)

    participantName = side
    meshName = side + "-Mesh"

    interface = precice.Interface(
        participantName, "../precice-config.xml", 0, 1)

    meshID = interface.get_mesh_id(meshName)

    writeData = "Temperature" if side == "Neumann" else "Flux"
    readData = "Flux" if side == "Neumann" else "Temperature"

    writedataID = interface.get_data_id(writeData, meshID)
    readdataID = interface.get_data_id(readData, meshID)

    couplinginterface = domain.boundary['right' if side ==
                                        'Dirichlet' else 'left']
    couplingsample = couplinginterface.sample('gauss', degree=degree * 2)

    vertices = couplingsample.eval(ns.x)

    dataIndices = interface.set_mesh_vertices(meshID, vertices)

    precice_dt = interface.initialize()

    if interface.is_action_required(action_write_initial_data()):
        writedata = numpy.zeros(len(dataIndices))
        interface.write_block_scalar_data(writedataID, dataIndices, writedata)
        interface.mark_action_fulfilled(action_write_initial_data())

    interface.initialize_data()

    t = 0

    projection = domain.integral(
        '(u - (1 + x_0 x_0 + alpha x_1 x_1 + beta ?t))^2' @ ns, degree=degree*2)
    lhs0 = nutils.solver.optimize(
        'lhs', projection, droptol=1e-15, arguments=dict(t=t))
    bezier = domain.sample('bezier', degree * 2)
    x, u, uexact = bezier.eval(['x_i', 'u', 'uexact'] @ ns, lhs=lhs0, t=t)

    with treelog.add(treelog.DataLog()):
        nutils.export.vtk(side + '-0',
                          bezier.tri, x, Temperature=u, reference=uexact)
                          
    t += precice_dt
    timestep = 1

    cons0 = nutils.solver.optimize(
        'lhs', sqr0, droptol=1e-15, arguments=dict(t=t))

    projectionmatrix = couplinginterface.integrate(
        nutils.function.outer(ns.basis), degree=degree * 2)
    projectioncons = numpy.zeros(res0.shape)
    projectioncons[projectionmatrix.rowsupp(1e-15)] = numpy.nan

    def fluxdofs(v): return projectionmatrix.solve(
        v, constrain=projectioncons)

    while interface.is_coupling_ongoing():
        if interface.is_action_required(action_write_iteration_checkpoint()):
            lhscheckpoint = lhs0
            conscheckpoint = cons0
            interface.mark_action_fulfilled(
                action_write_iteration_checkpoint())

        if interface.is_read_data_available():
            readdata = interface.read_block_scalar_data(
                readdataID, dataIndices)
            coupledata = couplingsample.asfunction(readdata)

            if side == 'Dirichlet':
                sqr = couplingsample.integral((ns.u - coupledata)**2)
                cons = nutils.solver.optimize(
                    'lhs', sqr, droptol=1e-15, constrain=cons0, arguments=dict(t=t))
                res = res0
            else:
                cons = cons0
                res = res0 + couplingsample.integral(ns.basis * coupledata / (nelems-1))

        dt = min(timestepsize, precice_dt)

        lhs = nutils.solver.solve_linear(
            'lhs', res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt, t=t))

        if interface.is_write_data_required(dt):
            if side == 'Dirichlet':
                fluxvalues = res.eval(lhs0=lhs0, lhs=lhs, dt=dt, t=t)
                writedata = couplingsample.eval(
                    'flux' @ ns, fluxdofs=fluxdofs(fluxvalues)) * (nelems-1)
            else:
                writedata = couplingsample.eval('u' @ ns, lhs=lhs)

            interface.write_block_scalar_data(
                writedataID, dataIndices, writedata)

        precice_dt = interface.advance(dt)

        if interface.is_action_required(action_read_iteration_checkpoint()):
            lhs0 = lhscheckpoint
            cons0 = conscheckpoint
            interface.mark_action_fulfilled(action_read_iteration_checkpoint())
        else:
            # prepare for next timestep
            lhs0 = lhs
            cons0 = nutils.solver.optimize(
                'lhs', sqr0, droptol=1e-15, arguments=dict(t=t+dt))

            if interface.is_time_window_complete():
                bezier = domain.sample('bezier', degree * 2)
                x, u, uexact = bezier.eval(
                    ['x_i', 'u', 'uexact'] @ ns, lhs=lhs0, t=t)

                with treelog.add(treelog.DataLog()):
                    nutils.export.vtk(side + "-" + str(timestep),
                                      bezier.tri, x, Temperature=u, reference=uexact)
                                      
            t += dt
            timestep += 1

    interface.finalize()


if __name__ == '__main__':
    nutils.cli.run(main)
