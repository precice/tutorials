import numpy as np
import matplotlib.pyplot as plt
import treelog
from nutils import mesh, function, solver, cli
import precice

def main(nelems=200, dt=.005, refdensity=1e3, refpressure=101325.0, psi=1e-6, viscosity=1e-3, theta=0.5):

    # --- preCICE initialization ---

    participant_name = "Fluid1DLeft"
    config_file_name = "../precice-config.xml"
    solver_process_index = 0
    solver_process_size = 1
    participant = precice.Participant(participant_name, config_file_name, solver_process_index, solver_process_size)
    mesh_name = "Fluid1DLeft-Mesh"
    velocity_name = "Velocity"
    pressure_name = "Pressure"
    positions = [[0.0, 500.0, 0.0]]     
    vertex_ids = participant.set_mesh_vertices(mesh_name, positions)

    participant.initialize()
    precice_dt = participant.get_max_time_step_size()

    # --- Nutils domain setup ---

    domain, geom = mesh.rectilinear([np.linspace(0, 500, nelems + 1)])    # 1D domain from 0 to 500 with nelems elements

    ns = function.Namespace()
    ns.x = geom
    ns.dt = dt
    ns.θ = theta
    ns.ρref = refdensity
    ns.pref = refpressure
    ns.pin = 98100  # Inlet pressure (Pa)
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

    # Weakly enforced inlet pressure boundary condition
    res += domain.boundary['left'].integral(
        'pin ubasis_ni n_i d:x' @ ns,
        degree=4
    )

    # Initial condition
    lhs0 = np.zeros(res.shape)

    # Time-stepping
    t = 0.0
    timestep = 0
    bezier = domain.sample('bezier', 2)

    f = open("watchpoint.txt", "w")

    # --- Time loop with preCICE coupling ---
    while participant.is_coupling_ongoing():

        # Read velocity data from other solver via preCICE
        u_read = participant.read_data(mesh_name, velocity_name, vertex_ids, precice_dt)

        # Constrain outlet velocity to match coupled value
        stringintegral = f'(u_0 - ({u_read[0][1]}))^2 d:x'
        sqr = domain.boundary['right'].integral(stringintegral @ ns, degree=4)
        cons = solver.optimize('lhs', sqr, droptol=1e-14)

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
            x, p, u, ρ = bezier.eval(['x_i','p', 'u_i', 'ρ'] @ ns, arguments=dict(lhs=lhs))

        # Send pressure at the right boundary to the other solver
        write_press = [[p[-1]]]
        participant.write_data(mesh_name, pressure_name, vertex_ids, write_press)

        # Advance in pseudo-time
        participant.advance(precice_dt)
        precice_dt = participant.get_max_time_step_size()

        # Restore checkpoint if needed
        if participant.requires_reading_checkpoint():
            lhs0 = lhscheckpoint
        else:
            # Update old solution
            lhs0 = lhs
            timestep += timestep

            # Save probe values (time, inlet pressure, inlet velocity, outlet pressure, outlet velocity, pressure at the middle, velocity at the middle)
            x, p, ρ, u = bezier.eval(['x_i', 'p', 'ρ', 'u_i'] @ ns, lhs=lhs)
            f.write("%e; %e; %e; %e; %e; %e; %e\n" % (t, p[0], u[0], p[-1], u[-1], p[199], u[199]))
            f.flush()

            t += precice_dt  # advance pseudo-time

    # Finalize preCICE
    participant.finalize()
    f.close()

if __name__ == '__main__':
    cli.run(main)