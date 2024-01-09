#! /usr/bin/env python3

from nutils import mesh, function, solver, export, cli
import numpy
import treelog
import precice

# for details on this solver see https://doi.org/10.1002/nme.6443

def main(inflow=10., viscosity=1., density=1., theta=.5, timestepsize=.01, npoints_per_elem=3):

    # mesh and geometry definition
    topo, geom = mesh.rectilinear([
        numpy.concatenate((
            numpy.linspace(-3, -1, 6, endpoint=False),
            numpy.linspace(-1, -0.3, 7, endpoint=False),
            numpy.linspace(-0.3, 0.3, 12, endpoint=False),
            numpy.linspace(0.3, 1, 7, endpoint=False),
            numpy.linspace(1, 3, 7))),
        numpy.concatenate((
            numpy.linspace(0, 1.5, 15, endpoint=False),
            numpy.linspace(1.5, 2, 3, endpoint=False),
            numpy.linspace(2, 4, 7))),
    ])

    domain = topo.withboundary(inflow='left', wall='top,bottom', outflow='right') - \
        topo[18:20, :10].withboundary(flap='left,right,top')

    couplinginterface = domain.boundary['flap']
    couplingsample = couplinginterface.sample('uniform', degree=npoints_per_elem)

    bezier = domain.sample('bezier', 2)

    # time approximations
    t0 = lambda f: function.replace_arguments(f, {arg: function.Argument(arg+'0', shape=shape, dtype=dtype)
                        for arg, (shape, dtype) in f.arguments.items() if arg in ('lhs', 'meshdofs', 'F')})
    # TR interpolation
    tθ = lambda f: theta * f + (1 - theta) * t0(f)
    # 1st order FD
    δt = lambda f: (f - t0(f)) / function.Argument('dt', ())

    # Nutils namespace
    ns = function.Namespace()
    ns.nu = viscosity
    ns.rho = density
    ns.uin = inflow
    ns.x0 = geom  # reference geometry
    ns.dbasis = domain.basis('std', degree=1).vector(2)
    ns.d_i = 'dbasis_ni ?meshdofs_n'
    ns.umesh_i = 'dbasis_ni ?umesh_n' # mesh velocity
    ns.x_i = 'x0_i + d_i'  # moving geometry
    ns.ubasis, ns.pbasis = function.chain([domain.basis('std', degree=2).vector(2), domain.basis('std', degree=1), ])
    ns.F_i = 'ubasis_ni ?F_n'  # stress field
    ns.urel_i = 'ubasis_ni ?lhs_n'  # relative velocity
    ns.u_i = 'umesh_i + urel_i'  # total velocity
    ns.p = 'pbasis_n ?lhs_n'  # pressure
    ns.qw = 1 / npoints_per_elem

    # boundary conditions for fluid equations
    sqr = domain.boundary['wall,flap'].integral('urel_k urel_k d:x0' @ ns, degree=4)
    cons = solver.optimize('lhs', sqr, droptol=1e-15)
    sqr = domain.boundary['inflow'].integral('((urel_0 - uin)^2 + urel_1^2) d:x0' @ ns, degree=4)
    cons = solver.optimize('lhs', sqr, droptol=1e-15, constrain=cons)

    # weak form fluid equations
    res = tθ(domain.integral('ubasis_ni,j (u_i,j + u_j,i) rho nu d:x' @ ns, degree=4))
    res += domain.integral('(-ubasis_ni,j p δ_ij + pbasis_n u_k,k) d:x' @ ns, degree=4)
    res += δt(domain.integral('rho ubasis_ni u_i d:x' @ ns, degree=4))
    res += tθ(domain.integral('rho ubasis_ni u_i,j urel_j d:x' @ ns, degree=4))

    # weak form for force computation
    resF = res + couplinginterface.sample('gauss', 4).integral('ubasis_ni F_i d:x' @ ns)
    consF = couplingsample.integrate((ns.ubasis**2).sum(1)) == 0

    # boundary conditions mesh displacements
    sqr = domain.boundary['inflow,outflow,wall'].integral('d_i d_i' @ ns, degree=2)
    meshcons0 = solver.optimize('meshdofs', sqr, droptol=1e-15)
    sqr = couplingsample.integral('d_k d_k' @ ns)
    meshcons = solver.optimize('meshdofs', sqr, droptol=1e-15, constrain=meshcons0)

    # weak form mesh displacements
    meshsqr = domain.integral('d_i,x0_j d_i,x0_j d:x0' @ ns, degree=2)

    # better initial guess: start from Stokes solution, comment out for comparison with other solvers
    #res_stokes = domain.integral('(ubasis_ni,j ((u_i,j + u_j,i) rho nu - p δ_ij) + pbasis_n u_k,k) d:x' @ ns, degree=4)
    #lhs0 = solver.solve_linear('lhs', res_stokes, constrain=cons, arguments=dict(meshdofs=meshdofs, meshdofs0=meshdofs0, meshdofs00=meshdofs00, meshdofs000=meshdofs000, dt=dt))

    # preCICE setup
    solverProcessIndex = 0
    solverProcessSize = 1
    interface = precice.Interface("Fluid", "../precice-config.xml", solverProcessIndex, solverProcessSize)
    meshID = interface.get_mesh_id("Fluid-Mesh")
    dataIndices = interface.set_mesh_vertices(meshID, couplingsample.eval(ns.x0))
    writedataID = interface.get_data_id("Force", meshID)
    readdataID = interface.get_data_id("Displacement", meshID)

    # initialize preCICE
    precice_dt = interface.initialize()

    timestep = 0
    arguments = dict(lhs=numpy.zeros(len(ns.ubasis)), meshdofs=numpy.zeros(len(ns.dbasis)))

    while interface.is_coupling_ongoing():
      with treelog.context(f'timestep {timestep}'):

        # read displacements from interface
        if interface.is_read_data_available():
            readdata = interface.read_block_vector_data(readdataID, dataIndices)
            coupledata = couplingsample.asfunction(readdata)
            sqr = couplingsample.integral(((ns.d - coupledata)**2).sum(0))
            meshcons = solver.optimize('meshdofs', sqr, droptol=1e-15, constrain=meshcons0)

        # save checkpoint
        if interface.is_action_required(precice.action_write_iteration_checkpoint()):
            checkpoint = timestep, arguments
            interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

        # advance variables
        timestep += 1
        arguments = {k+'0': arguments[k] for k in ('lhs', 'meshdofs')}
        arguments['dt'] = dt = min(precice_dt, timestepsize)
        arguments['meshdofs'] = solver.optimize('meshdofs', meshsqr, constrain=meshcons) # solve mesh deformation
        arguments['umesh'] = (arguments['meshdofs'] - arguments['meshdofs0']) / dt
        arguments['lhs'] = arguments['lhs0'] # initial guess for newton
        arguments['lhs'] = solver.newton('lhs', res, arguments=arguments, constrain=cons).solve(tol=1e-6) # solve fluid equations
        arguments['F'] = solver.solve_linear('F', resF, constrain=consF, arguments=arguments)

        # write forces to interface
        if interface.is_write_data_required(dt):
            writedata = couplingsample.eval('F_i qw d:x' @ ns, **arguments)
            interface.write_block_vector_data(writedataID, dataIndices, writedata)

        # do the coupling
        precice_dt = interface.advance(dt)

        # read checkpoint if required
        if interface.is_action_required(precice.action_read_iteration_checkpoint()):
            timestep, arguments = checkpoint
            interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())

        if interface.is_time_window_complete():
            x, u, p = bezier.eval(['x_i', 'u_i', 'p'] @ ns, **arguments)
            export.triplot('velocity.jpg', x, numpy.linalg.norm(u, axis=1), tri=bezier.tri, cmap='jet')
            export.triplot('pressure.jpg', x, p, tri=bezier.tri, cmap='jet')
            with treelog.add(treelog.DataLog()):
                export.vtk('Fluid_' + str(timestep), bezier.tri, x, u=u, p=p)

    interface.finalize()


if __name__ == '__main__':
    cli.run(main)
