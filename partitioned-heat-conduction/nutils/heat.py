#! /usr/bin/env python3

import nutils
import numpy
import treelog
import precice


def main(nelems: 'number of elements along edge' = 10, 
         etype: 'type of elements (square/triangle/mixed)' = 'square',
         btype: 'type of basis function (std/spline)' = 'std', 
         degree: 'polynomial degree' = 1, 
         side='Dirichlet', 
         timestepsize=.1):

    print("Running utils")

    domain, geom = nutils.mesh.unitsquare(nelems, etype)

    if side == 'Dirichlet':
        domain = domain[:nelems // 2, :]
    elif side == 'Neumann':
        domain = domain[nelems // 2:, :]
    else:
        raise Exception('invalid side {!r}'.format(side))

    ns = nutils.function.Namespace()
    ns.diffusivity = 1
    ns.x = geom
    ns.basis = domain.basis(btype, degree=degree)
    ns.u = 'basis_n ?lhs_n'
    ns.dudt = 'basis_n (?lhs_n - ?lhs0_n) / ?dt'
    ns.flux = 'basis_n ?fluxdofs_n'
    ns.uexact = 'sin(x_0) cosh(x_1)'

    res = domain.integral('(basis_n dudt + diffusivity basis_n,i u_,i) d:x' @ ns, degree=degree * 2)
    res -= domain.boundary['right'].integral('basis_n cos(1) cosh(x_1) d:x' @ ns, degree=degree * 2)

    sqr = domain.boundary['top'].integral('(u - cosh(1) sin(x_0))^2 d:x' @ ns, degree=degree * 2)
    if side == 'Dirichlet':
        sqr += domain.boundary['left'].integral('u^2 d:x' @ ns, degree=degree * 2)
        
    cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15)

    participantName = side
    meshName = side + "-Mesh"

    interface = precice.Interface(participantName, "../precice-config.xml", 0, 1)

    meshID = interface.get_mesh_id(meshName)

    writeData = "Temperature" if side == "Neumann" else "Flux"
    readData = "Flux" if side == "Neumann" else "Temperature"

    writedataID = interface.get_data_id(writeData, meshID)
    readdataID = interface.get_data_id(readData, meshID)

    couplinginterface = domain.boundary['right' if side == 'Dirichlet' else 'left']
    couplingsample = couplinginterface.sample('gauss', degree=degree * 2)

    vertices = couplingsample.eval(ns.x)

    dataIndices = interface.set_mesh_vertices(meshID, vertices)

    precice_dt = interface.initialize()
    
    if interface.is_action_required(precice.action_write_initial_data()):
      writedata = numpy.zeros(len(dataIndices))
      interface.write_block_scalar_data(writedataID, dataIndices, writedata)
      interface.mark_action_fulfilled(precice.action_write_initial_data())
    
    interface.initialize_data()

    cons0 = cons
    res0 = res

    projectionmatrix = couplinginterface.integrate(nutils.function.outer(ns.basis), degree=degree * 2)
    projectioncons = numpy.zeros(res.shape)
    projectioncons[projectionmatrix.rowsupp(1e-15)] = numpy.nan
    def fluxdofs(v): return projectionmatrix.solve(v, constrain=projectioncons)

    lhs0 = numpy.zeros(res.shape)

    while interface.is_coupling_ongoing():

        if interface.is_read_data_available():
            readdata = interface.read_block_scalar_data(readdataID, dataIndices)
            coupledata = couplingsample.asfunction(readdata)

            if side == 'Dirichlet':
                sqr = couplingsample.integral((ns.u - coupledata)**2)
                cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15, constrain=cons0)
            else:
                res = res0 + couplingsample.integral(ns.basis * coupledata)

        if interface.is_action_required(precice.action_write_iteration_checkpoint()):
            lhscheckpoint = lhs0
            interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())
            

        dt = min(timestepsize, precice_dt)

        lhs = nutils.solver.solve_linear(
            'lhs', res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt))

        if interface.is_write_data_required(dt):
            if side == 'Dirichlet':
                fluxvalues = res.eval(lhs0=lhs0, lhs=lhs, dt=dt)
                writedata = couplingsample.eval('flux' @ ns, fluxdofs=fluxdofs(fluxvalues))
            else:
                writedata = couplingsample.eval('u' @ ns, lhs=lhs)

            interface.write_block_scalar_data(writedataID, dataIndices, writedata)

        precice_dt = interface.advance(dt)
        
        lhs0 = lhs

        if interface.is_action_required(precice.action_read_iteration_checkpoint()):
            lhs0 = lhscheckpoint
            interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())
            
        if interface.is_time_window_complete():    
            bezier = domain.sample('bezier', degree * 2)
            x, u, uexact = bezier.eval(['x_i', 'u', 'uexact'] @ ns, lhs=lhs0)
            with treelog.add(treelog.DataLog()):
                nutils.export.vtk('solution-' + side, bezier.tri, x, fem=u, exact=uexact)

    interface.finalize()

if __name__ == '__main__':
    nutils.cli.run(main)
