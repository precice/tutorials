#! /usr/bin/env python3

from nutils import mesh, function, solver, export, cli
import numpy
import treelog
import precice

# for details on this solver see https://doi.org/10.1002/nme.6443

# some helper function to shift variables by one timestep


@function.replace
def subs0(f):
    if isinstance(f, function.Argument) and f._name == 'lhs':
        return function.Argument(name='lhs0', shape=f.shape, nderiv=f._nderiv)
    if isinstance(f, function.Argument) and f._name == 'meshdofs':
        return function.Argument(name='oldmeshdofs', shape=f.shape, nderiv=f._nderiv)
    if isinstance(f, function.Argument) and f._name == 'oldmeshdofs':
        return function.Argument(name='oldoldmeshdofs', shape=f.shape, nderiv=f._nderiv)
    if isinstance(f, function.Argument) and f._name == 'oldoldmeshdofs':
        return function.Argument(name='oldoldoldmeshdofs', shape=f.shape, nderiv=f._nderiv)

# some helper function to shift variables by two timesteps


@function.replace
def subs00(f):
    if isinstance(f, function.Argument) and f._name == 'lhs':
        return function.Argument(name='lhs00', shape=f.shape, nderiv=f._nderiv)


def main(inflow: 'inflow velocity' = 10,
         viscosity: 'kinematic viscosity' = 1.0,
         density: 'density' = 1.0,
         theta=0.5,
         timestepsize=0.01):

    # mesh and geometry definition
    grid_x_1 = numpy.linspace(-3, -1, 7)
    grid_x_1 = grid_x_1[:-1]
    grid_x_2 = numpy.linspace(-1, -0.3, 8)
    grid_x_2 = grid_x_2[:-1]
    grid_x_3 = numpy.linspace(-0.3, 0.3, 13)
    grid_x_3 = grid_x_3[:-1]
    grid_x_4 = numpy.linspace(0.3, 1, 8)
    grid_x_4 = grid_x_4[:-1]
    grid_x_5 = numpy.linspace(1, 3, 7)
    grid_x = numpy.concatenate((grid_x_1, grid_x_2, grid_x_3, grid_x_4, grid_x_5), axis=None)
    grid_y_1 = numpy.linspace(0, 1.5, 16)
    grid_y_1 = grid_y_1[:-1]
    grid_y_2 = numpy.linspace(1.5, 2, 4)
    grid_y_2 = grid_y_2[:-1]
    grid_y_3 = numpy.linspace(2, 4, 7)
    grid_y = numpy.concatenate((grid_y_1, grid_y_2, grid_y_3), axis=None)
    grid = [grid_x, grid_y]

    topo, geom = mesh.rectilinear(grid)
    domain = topo.withboundary(inflow='left', wall='top,bottom', outflow='right') - \
        topo[18:20, :10].withboundary(flap='left,right,top')

    # Nutils namespace
    ns = function.Namespace()

    # time approximations
    # TR interpolation
    ns._functions['t'] = lambda f: theta * f + (1 - theta) * subs0(f)
    ns._functions_nargs['t'] = 1
    # 1st order FD
    ns._functions['δt'] = lambda f: (f - subs0(f)) / dt
    ns._functions_nargs['δt'] = 1
    # 2nd order FD
    ns._functions['tt'] = lambda f: (1.5 * f - 2 * subs0(f) + 0.5 * subs00(f)) / dt
    ns._functions_nargs['tt'] = 1
    # extrapolation for pressure
    ns._functions['tp'] = lambda f: (1.5 * f - 0.5 * subs0(f))
    ns._functions_nargs['tp'] = 1

    ns.nu = viscosity
    ns.rho = density
    ns.uin = inflow
    ns.x0 = geom  # reference geometry
    ns.dbasis = domain.basis('std', degree=1).vector(2)
    ns.d_i = 'dbasis_ni ?meshdofs_n'
    ns.umesh_i = 'dbasis_ni (1.5 ?meshdofs_n - 2 ?oldmeshdofs_n + 0.5 ?oldoldmeshdofs_n ) / ?dt'
    ns.x_i = 'x0_i + d_i'  # moving geometry
    ns.ubasis, ns.pbasis = function.chain([domain.basis('std', degree=2).vector(2), domain.basis('std', degree=1), ])
    ns.F_i = 'ubasis_ni ?F_n'  # stress field
    ns.urel_i = 'ubasis_ni ?lhs_n'  # relative velocity
    ns.u_i = 'umesh_i + urel_i'  # total velocity
    ns.p = 'pbasis_n ?lhs_n'  # pressure

    # initialization of dofs
    meshdofs = numpy.zeros(len(ns.dbasis))
    oldmeshdofs = meshdofs
    oldoldmeshdofs = meshdofs
    oldoldoldmeshdofs = meshdofs
    lhs0 = numpy.zeros(len(ns.ubasis))

    # for visualization
    bezier = domain.sample('bezier', 2)

    # preCICE setup
    configFileName = "../precice-config.xml"
    participantName = "Fluid"
    solverProcessIndex = 0
    solverProcessSize = 1
    interface = precice.Interface(participantName, configFileName, solverProcessIndex, solverProcessSize)

    # define coupling meshes
    meshName = "Fluid-Mesh"
    meshID = interface.get_mesh_id(meshName)

    couplinginterface = domain.boundary['flap']
    couplingsample = couplinginterface.sample('gauss', degree=2)  # mesh located at Gauss points
    dataIndices = interface.set_mesh_vertices(meshID, couplingsample.eval(ns.x0))

    # coupling data
    writeData = "Force"
    readData = "Displacement"
    writedataID = interface.get_data_id(writeData, meshID)
    readdataID = interface.get_data_id(readData, meshID)

    # initialize preCICE
    precice_dt = interface.initialize()
    dt = min(precice_dt, timestepsize)

    # boundary conditions for fluid equations
    sqr = domain.boundary['wall,flap'].integral('urel_k urel_k d:x0' @ ns, degree=4)
    cons = solver.optimize('lhs', sqr, droptol=1e-15)
    sqr = domain.boundary['inflow'].integral('((urel_0 - uin)^2 + urel_1^2) d:x0' @ ns, degree=4)
    cons = solver.optimize('lhs', sqr, droptol=1e-15, constrain=cons)

    # weak form fluid equations
    res = domain.integral('t(ubasis_ni,j (u_i,j + u_j,i) rho nu d:x)' @ ns, degree=4)
    res += domain.integral('(-ubasis_ni,j p δ_ij + pbasis_n u_k,k) d:x' @ ns, degree=4)
    res += domain.integral('rho ubasis_ni δt(u_i d:x)' @ ns, degree=4)
    res += domain.integral('rho ubasis_ni t(u_i,j urel_j d:x)' @ ns, degree=4)

    # weak form for force computation
    resF = domain.integral('(ubasis_ni,j (u_i,j + u_j,i) rho nu d:x)' @ ns, degree=4)
    resF += domain.integral('tp(-ubasis_ni,j p δ_ij d:x)' @ ns, degree=4)
    resF += domain.integral('pbasis_n u_k,k d:x' @ ns, degree=4)
    resF += domain.integral('rho ubasis_ni tt(u_i d:x)' @ ns, degree=4)
    resF += domain.integral('rho ubasis_ni (u_i,j urel_j d:x)' @ ns, degree=4)
    resF += couplinginterface.sample('gauss', 4).integral('ubasis_ni F_i d:x' @ ns)
    consF = numpy.isnan(solver.optimize('F', couplinginterface.sample('gauss', 4).integral('F_i F_i' @ ns),
                                        droptol=1e-10))

    # boundary conditions mesh displacements
    sqr = domain.boundary['inflow,outflow,wall'].integral('d_i d_i' @ ns, degree=2)
    meshcons0 = solver.optimize('meshdofs', sqr, droptol=1e-15)

    # weak form mesh displacements
    meshsqr = domain.integral('d_i,x0_j d_i,x0_j d:x0' @ ns, degree=2)

    # better initial guess: start from Stokes solution, comment out for comparison with other solvers
    #res_stokes = domain.integral('(ubasis_ni,j ((u_i,j + u_j,i) rho nu - p δ_ij) + pbasis_n u_k,k) d:x' @ ns, degree=4)
    #lhs0 = solver.solve_linear('lhs', res_stokes, constrain=cons, arguments=dict(meshdofs=meshdofs, oldmeshdofs=oldmeshdofs, oldoldmeshdofs=oldoldmeshdofs, oldoldoldmeshdofs=oldoldoldmeshdofs, dt=dt))
    lhs00 = lhs0

    timestep = 0
    t = 0

    while interface.is_coupling_ongoing():

        # read displacements from interface
        if interface.is_read_data_available():
            readdata = interface.read_block_vector_data(readdataID, dataIndices)
            coupledata = couplingsample.asfunction(readdata)
            sqr = couplingsample.integral(((ns.d - coupledata)**2).sum(0))
            meshcons = solver.optimize('meshdofs', sqr, droptol=1e-15, constrain=meshcons0)
            meshdofs = solver.optimize('meshdofs', meshsqr, constrain=meshcons)

        # save checkpoint
        if interface.is_action_required(precice.action_write_iteration_checkpoint()):
            lhs_checkpoint = lhs0
            lhs00_checkpoint = lhs00
            t_checkpoint = t
            timestep_checkpoint = timestep
            oldmeshdofs_checkpoint = oldmeshdofs
            oldoldmeshdofs_checkpoint = oldoldmeshdofs
            oldoldoldmeshdofs_checkpoint = oldoldoldmeshdofs
            interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

        # solve fluid equations
        lhs1 = solver.newton('lhs', res, lhs0=lhs0, constrain=cons,
                             arguments=dict(lhs0=lhs0, dt=dt, meshdofs=meshdofs, oldmeshdofs=oldmeshdofs,
                                            oldoldmeshdofs=oldoldmeshdofs, oldoldoldmeshdofs=oldoldoldmeshdofs)
                             ).solve(tol=1e-6)

        # write forces to interface
        if interface.is_write_data_required(dt):
            F = solver.solve_linear('F', resF, constrain=consF,
                                    arguments=dict(lhs00=lhs00, lhs0=lhs0, lhs=lhs1, dt=dt, meshdofs=meshdofs,
                                                   oldmeshdofs=oldmeshdofs, oldoldmeshdofs=oldoldmeshdofs,
                                                   oldoldoldmeshdofs=oldoldoldmeshdofs))
            # writedata = couplingsample.eval(ns.F, F=F) # for stresses
            writedata = couplingsample.eval('F_i d:x' @ ns, F=F, meshdofs=meshdofs) * \
                numpy.concatenate([p.weights for p in couplingsample.points])[:, numpy.newaxis]
            interface.write_block_vector_data(writedataID, dataIndices, writedata)

        # do the coupling
        precice_dt = interface.advance(dt)
        dt = min(precice_dt, timestepsize)

        # advance variables
        timestep += 1
        t += dt
        lhs00 = lhs0
        lhs0 = lhs1
        oldoldoldmeshdofs = oldoldmeshdofs
        oldoldmeshdofs = oldmeshdofs
        oldmeshdofs = meshdofs

        # read checkpoint if required
        if interface.is_action_required(precice.action_read_iteration_checkpoint()):
            lhs0 = lhs_checkpoint
            lhs00 = lhs00_checkpoint
            t = t_checkpoint
            timestep = timestep_checkpoint
            oldmeshdofs = oldmeshdofs_checkpoint
            oldoldmeshdofs = oldoldmeshdofs_checkpoint
            oldoldoldmeshdofs = oldoldoldmeshdofs_checkpoint
            interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())

        if interface.is_time_window_complete():
            x, u, p = bezier.eval(['x_i', 'u_i', 'p'] @ ns, lhs=lhs1, meshdofs=meshdofs, oldmeshdofs=oldmeshdofs,
                                  oldoldmeshdofs=oldoldmeshdofs, oldoldoldmeshdofs=oldoldoldmeshdofs, dt=dt)
            with treelog.add(treelog.DataLog()):
                export.vtk('Fluid_' + str(timestep), bezier.tri, x, u=u, p=p)

    interface.finalize()


if __name__ == '__main__':
    cli.run(main)
