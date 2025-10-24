import numpy as np
import matplotlib.pyplot as plt
import treelog
import precice
from nutils import mesh, function, solver, cli

def main(nelems=200, dt=0.01, theta=0.5, u_in=0.1, rho=1, nu=10, R=5):
    """
    1D fluid solver (left domain) coupled via preCICE.
    Sends velocity to the right solver and reads pressure from it.
    """

    # --- preCICE setup ---
    participant_name = "Fluid1DLeft"
    config_file_name = "../precice-config.xml"
    participant = precice.Participant(participant_name, config_file_name, 0, 1)

    mesh_name = "Fluid1DLeft-Mesh"
    velocity_name = "Velocity"
    pressure_name = "Pressure"

    # Define coupling point (interface)
    positions = [[0.0, 0.0, 20.0]]
    vertex_ids = participant.set_mesh_vertices(mesh_name, positions)
    participant.initialize()
    precice_dt = participant.get_max_time_step_size()

    # --- Build 1D domain (from 0 to 20) ---
    domain, geom = mesh.rectilinear([np.linspace(0, 20, nelems + 1)])

    # --- Define parameters and variables ---
    ns = function.Namespace()
    ns.x = geom
    ns.dt = dt
    ns.θ = theta
    ns.ρ = rho
    ns.ν = nu
    ns.R = R
    ns.uin = u_in
    ns.α = '8 ν / (R^2)'  # Simple wall friction term

    # Basis functions for velocity (quadratic) and pressure (linear)
    ns.ubasis, ns.pbasis = function.chain([
        domain.basis('std', degree=2).vector(domain.ndims),
        domain.basis('std', degree=1)
    ])

    # Define fields and θ-scheme
    ns.u_i = 'ubasis_ni ?lhs_n'
    ns.u0_i = 'ubasis_ni ?lhs0_n'
    ns.p = 'pbasis_n ?lhs_n'
    ns.p0 = 'pbasis_n ?lhs0_n'
    ns.utheta_i = 'θ u_i + (1 - θ) u0_i'

    # --- Governing equations (momentum + continuity) ---
    res = domain.integral(
        'ubasis_ni (ρ (u_i - u0_i) / dt + ρ utheta_i utheta_j,j + α utheta_i + p_,i) d:x' @ ns,
        degree=4
    )
    res += domain.integral('pbasis_n utheta_i,i d:x' @ ns, degree=4)

    # Boundary condition: fixed inlet velocity
    bc = domain.boundary['left'].integral('(u_0 - uin)^2 d:x' @ ns, degree=4)
    bc0 = bc

    lhs0 = np.zeros(res.shape)  # Initial state
    t = 0.0
    bezier = domain.sample('bezier', 2)
    f = open("watchpoint.txt", "w")

    # --- Time stepping and coupling loop ---
    while participant.is_coupling_ongoing():
        # Read pressure from the right solver
        p_read = participant.read_data(mesh_name, pressure_name, vertex_ids, precice_dt)

        # Apply outlet pressure (coupled boundary)
        stringintegral = f'(p - ({p_read[0]}))^2 d:x'
        bc = bc0 + domain.boundary['right'].integral(stringintegral @ ns, degree=4)
        cons = solver.optimize('lhs', bc, droptol=1e-14)

        # Store checkpoint if preCICE requests
        if participant.requires_writing_checkpoint():
            lhscheckpoint = lhs0

        # Solve nonlinear system
        with treelog.context('newton'):
            lhs = solver.newton('lhs', residual=res, constrain=cons,
                                arguments=dict(lhs0=lhs0), lhs0=lhs0).solve(tol=1e-8)
            x, p, u = bezier.eval(['x_i', 'p', 'u_i'] @ ns, lhs=lhs)

        # Write velocity at the interface back to preCICE
        write_vel = [[0, 0, u[-1][0]]]
        participant.write_data(mesh_name, velocity_name, vertex_ids, write_vel)
        participant.advance(precice_dt)
        precice_dt = participant.get_max_time_step_size()

        # Handle checkpoint restore if needed
        if participant.requires_reading_checkpoint():
            lhs0 = lhscheckpoint
        else:
            lhs0 = lhs
            # Log pressure and velocity values at a few positions
            f.write("%e; %e; %e; %e; %e; %e; %e\n" %
                    (t, p[0], u[0], p[-1], u[-1], p[199], u[199]))
            f.flush()
            t += precice_dt

    # --- Finalization ---
    participant.finalize()
    np.savetxt("final_fields.txt", np.column_stack([x, u, p]), header="x u p")
    f.close()


if __name__ == '__main__':
    cli.run(main)
