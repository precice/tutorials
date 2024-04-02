#! /usr/bin/env python3

from nutils import cli, mesh, function, solver, export
import functools
import treelog
import numpy as np
import precice


def main(n=10, degree=1, timestep=.1, alpha=3., beta=1.2):

    x_grid = np.linspace(0, 1, n)
    y_grid = np.linspace(0, 1, n)

    # define the Nutils mesh
    domain, geom = mesh.rectilinear([x_grid, y_grid])
    coupling_boundary = domain.boundary['right']
    coupling_sample = coupling_boundary.sample('gauss', degree=degree * 2)

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
    ns.readbasis = coupling_sample.basis()
    ns.readfunc = 'readbasis_n ?readdata_n'

    # define the weak form
    res = domain.integral(
        '(basis_n dudt - basis_n f + basis_n,i u_,i) d:x' @ ns,
        degree=degree * 2)

    # set boundary conditions at non-coupling boundaries
    # top and bottom boundary are non-coupling for both sides
    sqr = domain.boundary['top,bottom,left'].integral('(u - uexact)^2 d:x' @ ns, degree=degree * 2)

    sqr += coupling_sample.integral('(u - readfunc)^2 d:x' @ ns)

    # preCICE setup
    participant = precice.Participant('Dirichlet', "../precice-config.xml", 0, 1)
    mesh_name = "Dirichlet-Mesh"
    vertex_ids = participant.set_mesh_vertices(
        mesh_name, coupling_sample.eval(ns.x))
    precice_write = functools.partial(
        participant.write_data,
        mesh_name,
        "Heat-Flux",
        vertex_ids)
    precice_read = functools.partial(
        participant.read_data,
        mesh_name,
        "Temperature",
        vertex_ids)

    # helper functions to project heat flux to coupling boundary
    # To communicate the flux to the Neumann side we should not simply
    # evaluate u_,i n_i as this is an unbounded term leading to suboptimal
    # convergence. Instead we project ∀ v: ∫_Γ v flux = ∫_Γ v u_,i n_i and
    # evaluate flux. While the right-hand-side contains the same unbounded
    # term, we can use the strong identity du/dt - u_,ii = f to rewrite it
    # to ∫_Ω [v (du/dt - f) + v_,i u_,i] - ∫_∂Ω\Γ v u_,k n_k, in which we
    # recognize the residual and an integral over the exterior boundary.
    # While the latter still contains the problematic unbounded term, we
    # can use the fact that the flux is a known value at the top and bottom
    # via the Dirichlet boundary condition, and impose it as constraints.
    rightsqr = domain.boundary['right'].integral(
        'flux^2 d:x' @ ns, degree=degree * 2)
    rightcons = solver.optimize('fluxdofs', rightsqr, droptol=1e-10)
    # rightcons is NaN in dofs that are NOT supported on the right boundary
    fluxsqr = domain.boundary['right'].boundary['top,bottom'].integral(
        '(flux - uexact_,0)^2 d:x' @ ns, degree=degree * 2)
    fluxcons = solver.optimize('fluxdofs',
                               fluxsqr,
                               droptol=1e-10,
                               constrain=np.choose(np.isnan(rightcons),
                                                   [np.nan,
                                                    0.]))
    # fluxcons is NaN in dofs that are supported on ONLY the right boundary
    fluxres = coupling_sample.integral('basis_n flux d:x' @ ns) - res

    # write initial data
    if participant.requires_initial_data():
        precice_write(coupling_sample.eval(0.))

    participant.initialize()
    precice_dt = participant.get_max_time_step_size()
    solver_dt = timestep
    dt = min(precice_dt, solver_dt)

    t = 0.
    istep = 0

    # initial condition
    sqr0 = domain.integral('(u - uexact)^2' @ ns, degree=degree * 2)
    lhs = solver.optimize('lhs', sqr0, arguments=dict(t=t))
    bezier = domain.sample('bezier', degree * 2)

    while True:

        # generate output
        x, u, uexact = bezier.eval(['x_i', 'u', 'uexact'] @ ns, lhs=lhs, t=t)
        with treelog.add(treelog.DataLog()):
            export.vtk(
                "Dirichlet-" +
                str(istep),
                bezier.tri,
                x,
                Temperature=u,
                reference=uexact)

        if not participant.is_coupling_ongoing():
            break

        # save checkpoint
        if participant.requires_writing_checkpoint():
            checkpoint = lhs, t, istep

        # prepare next timestep
        precice_dt = participant.get_max_time_step_size()
        dt = min(solver_dt, precice_dt)
        lhs0 = lhs
        istep += 1
        # read data from participant
        readdata = precice_read(dt)
        t += dt

        # update (time-dependent) boundary condition
        cons = solver.optimize(
            'lhs',
            sqr,
            droptol=1e-15,
            arguments=dict(
                t=t,
                readdata=readdata))

        # solve nutils timestep
        lhs = solver.solve_linear(
            'lhs', res, constrain=cons, arguments=dict(
                lhs0=lhs0, dt=dt, t=t, readdata=readdata))

        # write data to participant
        fluxdofs = solver.solve_linear(
            'fluxdofs', fluxres, arguments=dict(
                lhs0=lhs0, lhs=lhs, dt=dt, t=t), constrain=fluxcons)
        write_data = coupling_sample.eval(
            'flux' @ ns, fluxdofs=fluxdofs)

        precice_write(write_data)

        # do the coupling
        participant.advance(dt)

        # read checkpoint if required
        if participant.requires_reading_checkpoint():
            lhs, t, istep = checkpoint

    participant.finalize()


if __name__ == '__main__':
    cli.run(main)
