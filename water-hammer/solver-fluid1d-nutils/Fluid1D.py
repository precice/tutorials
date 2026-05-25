import numpy as np
import matplotlib.pyplot as plt
import treelog
from nutils import mesh, function, solver, cli
import precice


def main(nelems=400, dt=.005, refdensity=1e3, refpressure=101325.0, psi=1e-6, viscosity=1e-3, theta=0.5, side='Left'):

    # --- preCICE initialization ---

    if side == "Left":
        participant_name = "Fluid1DLeft"
        mesh_name = "Fluid1DLeft-Mesh"
        read_data = "Velocity"
        write_data = "Pressure"
    elif side == "Right":
        participant_name = "Fluid1DRight"
        mesh_name = "Fluid1DRight-Mesh"
        read_data = "Pressure"
        write_data = "Velocity"
    else:
        raise Exception('invalid side {!r}'.format(side))
    config_file_name = "../precice-config.xml"
    solver_process_index = 0
    solver_process_size = 1
    participant = precice.Participant(participant_name, config_file_name, solver_process_index, solver_process_size)
    positions = [[0.0, 500.0, 0.0]]     # Define a single coupling point (3D format required by preCICE)
    vertex_ids = participant.set_mesh_vertices(mesh_name, positions)

    participant.initialize()
    precice_dt = participant.get_max_time_step_size()

    # --- Nutils domain setup ---

    domain, geom = mesh.rectilinear([np.linspace(0, 1000, nelems + 1)])

    if side == "Left":
        domain = domain[:nelems // 2]
    elif side == "Right":
        domain = domain[nelems // 2:]
    else:
        raise Exception('invalid side {!r}'.format(side))

    ns = function.Namespace()
    ns.x = geom
    ns.dt = dt
    ns.θ = theta
    ns.ρref = refdensity
    ns.pref = refpressure
    ns.pin = 98100
    ns.μ = viscosity
    ns.ψ = psi

    # Define basis functions: velocity (vector, degree 2), density (scalar, degree 1)
    ns.ubasis, ns.ρbasis = function.chain([
        domain.basis('std', degree=2).vector(domain.ndims),
        domain.basis('std', degree=1)
    ])

    # Define trial and previous-step functions
    ns.u_i = 'ubasis_ni ?lhs_n'
    ns.u0_i = 'ubasis_ni ?lhs0_n'
    ns.ρ = 'ρref + ρbasis_n ?lhs_n'
    ns.ρ0 = 'ρref + ρbasis_n ?lhs0_n'
    ns.p = 'pref + (ρ - ρref) / ψ'
    ns.utheta_i = 'θ u_i + (1 - θ) u0_i'
    ns.ρtheta = 'θ ρ + (1 - θ) ρ0'

    # Cauchy stress tensor: viscous + pressure
    ns.σ_ij = 'μ (utheta_i,j + utheta_j,i) - p δ_ij'

    # --- Weak form residuals ---

    # Mass conservation
    res = domain.integral(
        'ρbasis_n ((ρ - ρ0) / dt + (ρtheta utheta_k)_,k) d:x' @ ns,
        degree=4
    )

    # Momentum conservation: unsteady + convective + diffusive terms
    res += domain.integral(
        '(ubasis_ni (((ρ u_i - ρ0 u0_i) / dt) + (ρtheta utheta_i utheta_j)_,j) + ubasis_ni,j σ_ij) d:x' @ ns,
        degree=4
    )

    if side == "Left":    # Weakly enforced inlet pressure boundary condition
        res += domain.boundary['left'].integral('pin ubasis_ni n_i d:x' @ ns, degree=4)
    elif side == "Right":  # Outlet velocity boundary condition
        sqr = domain.boundary['right'].integral('(u_0 - 1)^2' @ ns, degree=4)
        cons0 = solver.optimize('lhs', sqr, droptol=1e-14)
        res0 = res        # Save base residual (without inlet BC)

    # Initial condition
    lhs0 = np.zeros(res.shape)

    # Time-stepping
    t = 0.0
    timestep = 0
    bezier = domain.sample('bezier', 2)

    f = open("probes.txt", "w")

    # --- Time loop with preCICE coupling ---
    while participant.is_coupling_ongoing():

        # Read data from other solver via preCICE
        data_read = participant.read_data(mesh_name, read_data, vertex_ids, precice_dt)

        if side == "Left":
            ns.uOut = data_read[0][1]
            sqr = domain.boundary['right'].integral('(u_0 - uOut)^2 d:x' @ ns, degree=4)
            cons = solver.optimize('lhs', sqr, droptol=1e-14)
        elif side == "Right":
            ns.pIn = data_read[0]
            res = res0 + domain.boundary['left'].integral('pIn ubasis_ni n_i d:x' @ ns, degree=4)
            cons = min(0.2 * t, 1) * cons0

        # Write checkpoint if required by preCICE
        if participant.requires_writing_checkpoint():
            lhscheckpoint = lhs0

        # Solve nonlinear system (Newton's method)
        with treelog.context('stokes'):
            lhs = solver.newton(
                'lhs',
                residual=res,
                constrain=cons,
                arguments=dict(lhs0=lhs0),
                lhs0=lhs0
            ).solve(1e-08)

            # Evaluate fields for visualization/debugging
            x, p, u, ρ = bezier.eval(['x_i', 'p', 'u_i', 'ρ'] @ ns, arguments=dict(lhs=lhs))

        # Send pressure at the right boundary to the other solver
        if side == "Left":
            value_written = [[p[-1]]]
        elif side == "Right":
            value_written = [[0, u[0][0], 0]]
        participant.write_data(mesh_name, write_data, vertex_ids, value_written)

        # Advance in pseudo-time
        participant.advance(precice_dt)
        precice_dt = participant.get_max_time_step_size()

        # Restore checkpoint if needed
        if participant.requires_reading_checkpoint():
            lhs0 = lhscheckpoint
        else:
            # Update old solution
            lhs0 = lhs
            timestep += 1

            # Save probe values (time, inlet pressure, inlet velocity, outlet
            # pressure, outlet velocity, pressure at the middle, velocity at the
            # middle)
            x, p, ρ, u = bezier.eval(['x_i', 'p', 'ρ', 'u_i'] @ ns, lhs=lhs)
            mid = len(p) // 2
            f.write("%e; %e; %e; %e; %e; %e; %e\n" % (t, p[0], u[0], p[-1], u[-1], p[mid], u[mid]))
            f.flush()

            t += precice_dt  # advance pseudo-time

    # Finalize preCICE
    participant.finalize()
    f.close()


if __name__ == '__main__':
    cli.run(main)
