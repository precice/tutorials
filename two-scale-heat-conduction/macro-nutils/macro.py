#! /usr/bin/env python3
#
# In this script we solve the unsteady Heat equation

from nutils import mesh, function, solver, export, cli
import treelog
import numpy as np
import precice


def main():
    """
    2D unsteady heat equation on a unit square.
    The material consists of a mixture of two materials, the grain and sand
    """
    is_coupled_case = True

    topo, geom = mesh.rectilinear([np.linspace(0, 1.0, 9), np.linspace(0, 0.5, 5)])

    ns = function.Namespace(fallback_length=2)
    ns.x = geom
    ns.basis = topo.basis('std', degree=1)
    ns.kbasis = topo.basis('std', degree=1).vector(topo.ndims).vector(topo.ndims)
    ns.u = 'basis_n ?solu_n'

    # Coupling quantities
    ns.phi = 'basis_n ?solphi_n'
    ns.k_ij = 'kbasis_nij ?solk_n'

    phi = 0.5  # initial value
    k = 1.0  # initial value

    ns.rhos = 1.0
    ns.rhog = 1.0
    ns.dudt = 'basis_n (?solu_n - ?solu0_n) / ?dt'

    # Dirichlet boundary condition on lower left corner
    ns.usource = 0.0

    # Initial condition
    ns.uinitial = 0.5

    if is_coupled_case:
        # preCICE setup
        interface = precice.Interface("Macro-heat", "precice-config-full-sim.xml", 0, 1)

        # define coupling meshes
        mesh_name = "macro-mesh"
        mesh_id = interface.get_mesh_id(mesh_name)

        # Define Gauss points on entire domain as coupling mesh
        couplingsample = topo.sample('gauss', degree=2)  # mesh located at Gauss points
        vertex_ids = interface.set_mesh_vertices(mesh_id, couplingsample.eval(ns.x))

        print("Number of coupling vertices = {}".format(len(vertex_ids)))

        # coupling data
        k_00_id = interface.get_data_id("k_00", mesh_id)
        k_01_id = interface.get_data_id("k_01", mesh_id)
        k_10_id = interface.get_data_id("k_10", mesh_id)
        k_11_id = interface.get_data_id("k_11", mesh_id)
        poro_id = interface.get_data_id("porosity", mesh_id)
        conc_id = interface.get_data_id("concentration", mesh_id)

        # initialize preCICE
        dt = interface.initialize()
    else:
        dt = 1.0E-2
        sqrphi = topo.integral((ns.phi - phi) ** 2, degree=1)
        solphi = solver.optimize('solphi', sqrphi, droptol=1E-12)

        sqrk = topo.integral(((ns.k - k * np.eye(2)) * (ns.k - k * np.eye(2))).sum([0, 1]), degree=2)
        solk = solver.optimize('solk', sqrk, droptol=1E-12)

    # Time related variables
    ns.dt = dt
    n = n_checkpoint = 0
    t = t_checkpoint = 0
    t_out = 0.05
    t_end = 0.25
    n_out = int(t_out / dt)
    n_t = int(t_end / dt)

    # define the weak form
    res = topo.integral('((rhos phi + (1 - phi) rhog) basis_n dudt + k_ij basis_n,i u_,j) d:x' @ ns, degree=2)

    # Set Dirichlet boundary conditions
    sqr = topo.boundary['left'].boundary['bottom'].integral('(u - usource)^2 d:x' @ ns, degree=2)
    cons = solver.optimize('solu', sqr, droptol=1e-15)

    # Set domain to initial condition
    sqr = topo.integral('(u - uinitial)^2' @ ns, degree=2)
    solu0 = solver.optimize('solu', sqr)

    if is_coupled_case:
        if interface.is_action_required(precice.action_write_initial_data()):
            temperatures = couplingsample.eval('u' @ ns, solu=solu0)
            interface.write_block_scalar_data(conc_id, vertex_ids, temperatures)

            interface.mark_action_fulfilled(precice.action_write_initial_data())

        interface.initialize_data()

    # Prepare the post processing sample
    bezier = topo.sample('bezier', 2)

    # VTK output of initial state
    x, u = bezier.eval(['x_i', 'u'] @ ns, solu=solu0)
    with treelog.add(treelog.DataLog()):
        export.vtk('macro-heat-0', bezier.tri, x, u=u)

    if is_coupled_case:
        is_coupling_ongoing = interface.is_coupling_ongoing()
    else:
        is_coupling_ongoing = True

    # time loop
    while is_coupling_ongoing:
        if is_coupled_case:
            # write checkpoint
            if interface.is_action_required(precice.action_write_iteration_checkpoint()):
                solu_checkpoint = solu0
                t_checkpoint = t
                n_checkpoint = n
                interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

            # Read porosity and apply
            poro_data = interface.read_block_scalar_data(poro_id, vertex_ids)
            poro_coupledata = couplingsample.asfunction(poro_data)
            sqrphi = couplingsample.integral((ns.phi - poro_coupledata) ** 2)
            solphi = solver.optimize('solphi', sqrphi, droptol=1E-12)

            # Read conductivity and apply
            k_00 = interface.read_block_scalar_data(k_00_id, vertex_ids)
            k_01 = interface.read_block_scalar_data(k_01_id, vertex_ids)
            k_10 = interface.read_block_scalar_data(k_10_id, vertex_ids)
            k_11 = interface.read_block_scalar_data(k_11_id, vertex_ids)

            k_00_c = couplingsample.asfunction(k_00)
            k_01_c = couplingsample.asfunction(k_01)
            k_10_c = couplingsample.asfunction(k_10)
            k_11_c = couplingsample.asfunction(k_11)

            k_coupledata = function.asarray([[k_00_c, k_01_c], [k_10_c, k_11_c]])
            sqrk = couplingsample.integral(((ns.k - k_coupledata) * (ns.k - k_coupledata)).sum([0, 1]))
            solk = solver.optimize('solk', sqrk, droptol=1E-12)

        # solve timestep
        solu = solver.solve_linear('solu', res, constrain=cons,
                                   arguments=dict(solu0=solu0, dt=dt, solphi=solphi, solk=solk))

        if is_coupled_case:
            concentrations = couplingsample.eval('u' @ ns, solu=solu)
            interface.write_block_scalar_data(conc_id, vertex_ids, concentrations)

            # do the coupling
            precice_dt = interface.advance(dt)
            dt = min(precice_dt, dt)

        # advance variables
        n += 1
        t += dt
        solu0 = solu

        if is_coupled_case:
            is_coupling_ongoing = interface.is_coupling_ongoing()

            if interface.is_action_required(precice.action_read_iteration_checkpoint()):
                solu0 = solu_checkpoint
                t = t_checkpoint
                n = n_checkpoint
                interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())
            else:
                if n % n_out == 0:
                    x, phi, k, u = bezier.eval(['x_i', 'phi', 'k', 'u'] @ ns, solphi=solphi, solk=solk, solu=solu)
                    with treelog.add(treelog.DataLog()):
                        export.vtk('macro-heat-' + str(n), bezier.tri, x, u=u, phi=phi, K=k)
        else:
            if n % n_out == 0:
                x, phi, k, u = bezier.eval(['x_i', 'phi', 'k', 'u'] @ ns, solphi=solphi, solk=solk, solu=solu)
                with treelog.add(treelog.DataLog()):
                    export.vtk('macro-heat-' + str(n), bezier.tri, x, u=u, phi=phi, K=k)

            if n >= n_t:
                is_coupling_ongoing = False

    if is_coupled_case:
        interface.finalize()


if __name__ == '__main__':
    cli.run(main)
