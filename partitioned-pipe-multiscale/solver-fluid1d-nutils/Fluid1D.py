import numpy as np
import matplotlib.pyplot as plt
import treelog
import precice
from nutils import mesh, function, solver, cli


def main(nelems=400, dt=0.01, theta=0.5, u_in=0.1, p_out=0, rho=1, nu=10, R=5, side="Left"):
    """
    1D fluid solver coupled via preCICE.
    """

    # --- preCICE setup ---

    if side == "Left":
        participant_name = "Fluid1DLeft"
        mesh_name = "Fluid1DLeft-Mesh"
        read_data = "Pressure"
        write_data = "Velocity"
    elif side == "Right":
        participant_name = "Fluid1DRight"
        mesh_name = "Fluid1DRight-Mesh"
        read_data = "Velocity"
        write_data = "Pressure"
    else:
        raise Exception('invalid side {!r}'.format(side))
    config_file_name = "../precice-config.xml"
    solver_process_index = 0
    solver_process_size = 1
    participant = precice.Participant(participant_name, config_file_name, solver_process_index, solver_process_size)
    positions = [[0.0, 0.0, 20.0]]  # Define a single coupling point
    vertex_ids = participant.set_mesh_vertices(mesh_name, positions)

    participant.initialize()
    precice_dt = participant.get_max_time_step_size()

    # --- Nutils domain setup ---

    domain, geom = mesh.rectilinear([np.linspace(0, 40, nelems + 1)])
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
    ns.ρ = rho
    ns.ν = nu
    ns.R = R
    ns.uin = u_in
    ns.pout = p_out
    ns.α = '8 ν / (R^2)'  # Simple wall friction term

    # Define basis functions: velocity (vector, degree 2), pressure (scalar, degree 1)
    ns.ubasis, ns.pbasis = function.chain([
        domain.basis('std', degree=2).vector(domain.ndims),
        domain.basis('std', degree=1)
    ])

    # Define trial and previous-step functions
    ns.u_i = 'ubasis_ni ?lhs_n'
    ns.u0_i = 'ubasis_ni ?lhs0_n'
    ns.p = 'pbasis_n ?lhs_n'
    ns.p0 = 'pbasis_n ?lhs0_n'
    ns.utheta_i = 'θ u_i + (1 - θ) u0_i'

    # --- Weak form residuals ---

    # Momentum conservation
    res = domain.integral(
        'ubasis_ni (ρ (u_i - u0_i) / dt + ρ utheta_i utheta_j,j + α utheta_i + p_,i) d:x' @ ns,
        degree=4
    )

    # Mass conservation
    res += domain.integral('pbasis_n utheta_i,i d:x' @ ns, degree=4)

    if side == "Left":
        bc0 = domain.boundary['left'].integral('(u_0 - uin)^2 d:x' @ ns, degree=4)
    elif side == "Right":
        bc0 = domain.boundary['right'].integral('(p - pout)^2 d:x' @ ns, degree=4)

    # Initial condition
    lhs0 = np.zeros(res.shape)

    # Time-stepping
    t = 0.0
    bezier = domain.sample('bezier', 2)

    f = open("probes.txt", "w")

    # --- Time loop with preCICE coupling ---
    while participant.is_coupling_ongoing():

        # Read data from other solver via preCICE
        data_read = participant.read_data(mesh_name, read_data, vertex_ids, precice_dt)

        if side == "Left":
            ns.pOut = data_read[0]
            bc = bc0 + domain.boundary['right'].integral('(p - pOut)^2 d:x' @ ns, degree=4)
        elif side == "Right":
            ns.uIn = data_read[0][2]
            bc = bc0 + domain.boundary['left'].integral('(u_0 - uIn)^2 d:x' @ ns, degree=4)
        cons = solver.optimize('lhs', bc, droptol=1e-14)

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
            ).solve(tol=1e-8)

        # Evaluate fields for visualization/debugging
        x, p, u = bezier.eval(['x_i', 'p', 'u_i'] @ ns, lhs=lhs)

        # Send data to 3D solver
        if side == "Left":
            value_written = [[0, 0, u[-1][0]]]
        elif side == "Right":
            value_written = [[p[0]]]
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
            # Save probe values (time, inlet pressure, inlet velocity, outlet
            # pressure, outlet velocity, pressure at the middle, velocity at the
            # middle)
            mid = len(p) // 2
            x, p, u = bezier.eval(['x_i', 'p', 'u_i'] @ ns, lhs=lhs)
            f.write("%e; %e; %e; %e; %e; %e; %e\n" %
                    (t, p[0], u[0], p[-1], u[-1], p[mid], u[mid]))

            f.flush()

            t += precice_dt  # Advance pseudo time

    # --- Finalization ---
    participant.finalize()
    np.savetxt("final_fields.txt", np.column_stack([x, u, p]), header="x u p")
    f.close()


if __name__ == '__main__':
    cli.run(main)
