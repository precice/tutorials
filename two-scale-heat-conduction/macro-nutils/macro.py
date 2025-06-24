#! /usr/bin/env python3
#
# Heat conduction problem with two materials.

from nutils import mesh, function, solver, export, cli
import treelog
import numpy as np
import precice


def main():
    """
    2D heat conduction problem solved on a rectangular domain.
    The material consists of a mixture of two materials "g" and "s".
    """
    is_coupled_case = True  # If False, single-physics problem is solved

    topo, geom = mesh.rectilinear([np.linspace(0, 1.0, 5), np.linspace(0, 0.5, 4)])

    ns = function.Namespace(fallback_length=2)
    ns.x = geom
    ns.basis = topo.basis('std', degree=1)
    ns.u = 'basis_n ?solu_n'  # the field to be solved for

    ns.phi = 'basis_n ?solphi_n'  # porosity, to be read from preCICE
    phi = 0.5  # initial value of porosity

    ns.kbasis = topo.basis('std', degree=1).vector(topo.ndims).vector(topo.ndims)
    # conductivity matrix, to be read from preCICE (read as two vectors and then patched to form a matrix)
    ns.k_ij = 'kbasis_nij ?solk_n'
    k = 1.0  # initial value of conductivity

    ns.rhos = 1.0  # density of material s
    ns.rhog = 1.0  # density of material g
    ns.dudt = 'basis_n (?solu_n - ?solu0_n) / ?dt'  # implicit Euler time stepping formulation

    # Dirichlet boundary condition at lower left corner
    ns.usource = 0.0

    # Initial condition in the domain
    ns.uinitial = 0.5

    if is_coupled_case:
        participant = precice.Participant("macro-heat", "../precice-config.xml", 0, 1)
        mesh_name = "macro-mesh"

        # Define Gauss points on entire domain as coupling mesh (volume coupling from macro side)
        couplingsample = topo.sample('gauss', degree=2)  # mesh vertices are Gauss points
        vertex_ids = participant.set_mesh_vertices(mesh_name, couplingsample.eval(ns.x))
    else:
        sqrphi = topo.integral((ns.phi - phi) ** 2, degree=1)
        solphi = solver.optimize('solphi', sqrphi, droptol=1E-12)

        sqrk = topo.integral(((ns.k - k * np.eye(2)) * (ns.k - k * np.eye(2))).sum([0, 1]), degree=2)
        solk = solver.optimize('solk', sqrk, droptol=1E-12)

    # Define the weak form of the heat conduction problem with two materials
    res = topo.integral('((rhos phi + (1 - phi) rhog) basis_n dudt + k_ij basis_n,i u_,j) d:x' @ ns, degree=2)

    # Set Dirichlet boundary conditions
    sqr = topo.boundary['left'].boundary['bottom'].integral('(u - usource)^2 d:x' @ ns, degree=2)
    cons = solver.optimize('solu', sqr, droptol=1e-15)

    # Set domain to initial condition
    sqr = topo.integral('(u - uinitial)^2' @ ns, degree=2)
    solu0 = solver.optimize('solu', sqr)

    if is_coupled_case:
        if participant.requires_initial_data():
            concentrations = couplingsample.eval('u' @ ns, solu=solu0)
            participant.write_data(mesh_name, "concentration", vertex_ids, concentrations)

        participant.initialize()
        dt = solver_dt = participant.get_max_time_step_size()
    else:
        dt = 1.0E-2

    ns.dt = dt
    n = n_checkpoint = 0
    t = t_checkpoint = 0
    t_out = 0.01
    t_end = 0.25  # Only relevant when single physics case is run
    n_out = int(t_out / dt)
    n_t = int(t_end / dt)

    # Prepare the post processing sample
    bezier = topo.sample('bezier', 2)

    # Output of initial state
    x, u = bezier.eval(['x_i', 'u'] @ ns, solu=solu0)
    with treelog.add(treelog.DataLog()):
        export.vtk('macro-0', bezier.tri, x, u=u)

    if is_coupled_case:
        is_coupling_ongoing = participant.is_coupling_ongoing()
    else:
        is_coupling_ongoing = True

    while is_coupling_ongoing:
        if is_coupled_case:
            if participant.requires_writing_checkpoint():
                solu_checkpoint = solu0
                t_checkpoint = t
                n_checkpoint = n

            precice_dt = participant.get_max_time_step_size()
            dt = min(precice_dt, solver_dt)

            # Read porosity and apply it to the existing solution
            poro_data = participant.read_data(mesh_name, "porosity", vertex_ids, dt)
            poro_coupledata = couplingsample.asfunction(poro_data)
            sqrphi = couplingsample.integral((ns.phi - poro_coupledata) ** 2)
            solphi = solver.optimize('solphi', sqrphi, droptol=1E-12)

            # Read conductivity and apply it to the existing solution
            k_00_c = couplingsample.asfunction(participant.read_data(mesh_name, "k_00", vertex_ids, dt))
            k_01_c = couplingsample.asfunction(participant.read_data(mesh_name, "k_01", vertex_ids, dt))
            k_10_c = couplingsample.asfunction(participant.read_data(mesh_name, "k_10", vertex_ids, dt))
            k_11_c = couplingsample.asfunction(participant.read_data(mesh_name, "k_11", vertex_ids, dt))

            conductivity = function.asarray([[k_00_c, k_01_c], [k_10_c, k_11_c]])
            sqrk = couplingsample.integral(((ns.k - conductivity) * (ns.k - conductivity)).sum([0, 1]))
            solk = solver.optimize('solk', sqrk, droptol=1E-12)

        solu = solver.solve_linear('solu', res, constrain=cons,
                                   arguments=dict(solu0=solu0, dt=dt, solphi=solphi, solk=solk))

        if is_coupled_case:
            # Collect values of field u and write them to preCICE
            concentration = couplingsample.eval('u' @ ns, solu=solu)
            participant.write_data(mesh_name, "concentration", vertex_ids, concentration)

            participant.advance(dt)

        n += 1
        t += dt
        solu0 = solu

        if is_coupled_case:
            is_coupling_ongoing = participant.is_coupling_ongoing()

            if participant.requires_reading_checkpoint():
                solu0 = solu_checkpoint
                t = t_checkpoint
                n = n_checkpoint
            else:
                if n % n_out == 0:
                    x, phi, k, u = bezier.eval(['x_i', 'phi', 'k', 'u'] @ ns, solphi=solphi, solk=solk, solu=solu)
                    with treelog.add(treelog.DataLog()):
                        export.vtk('macro-' + str(n), bezier.tri, x, u=u, phi=phi, K=k)
        else:
            if n % n_out == 0:
                x, phi, k, u = bezier.eval(['x_i', 'phi', 'k', 'u'] @ ns, solphi=solphi, solk=solk, solu=solu)
                with treelog.add(treelog.DataLog()):
                    export.vtk('macro-' + str(n), bezier.tri, x, u=u, phi=phi, K=k)

            if n >= n_t:
                is_coupling_ongoing = False

    if is_coupled_case:
        participant.finalize()


if __name__ == '__main__':
    cli.run(main)
