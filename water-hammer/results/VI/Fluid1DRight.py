import numpy as np
import matplotlib.pyplot as plt
import treelog
from nutils import mesh, function, solver, cli
import precice


def main(nelems=200, dt=.005, refdensity=1e3, refpressure=101325.0, psi=1e-6, viscosity=1e-3, theta=0.5):

    # Initialize preCICE
    participant_name = "Fluid1DRight"
    config_file_name = "../precice-config.xml"
    solver_process_index = 0
    solver_process_size = 1
    participant = precice.Participant(participant_name, config_file_name, solver_process_index, solver_process_size)
    mesh_name = "Fluid1DRight-Mesh"
    velocity_name = "Velocity"
    pressure_name = "Pressure"
    positions = [[0, 500.0, 0]]  # one vertex of three dimensions
    vertex_ids = participant.set_mesh_vertices(mesh_name, positions)
    participant.initialize()
    precice_dt = participant.get_max_time_step_size()

    domain, geom = mesh.rectilinear([np.linspace(500, 1000, nelems + 1)])

    ns = function.Namespace()
    ns.x = geom
    ns.dt = dt
    ns.θ = theta
    ns.ρref = refdensity
    ns.pref = refpressure
    ns.μ = viscosity
    ns.ψ = psi

    ns.ubasis, ns.ρbasis = function.chain(
        [domain.basis('std', degree=2).vector(domain.ndims), domain.basis('std', degree=1)])

    ns.u_i = 'ubasis_ni ?lhs_n'
    ns.u0_i = 'ubasis_ni ?lhs0_n'
    ns.ρ = 'ρref + ρbasis_n ?lhs_n'
    ns.ρ0 = 'ρref + ρbasis_n ?lhs0_n'
    ns.p = 'pref + (ρ - ρref) / ψ'
    ns.utheta_i = 'θ u_i + (1 - θ) u0_i'
    ns.ρtheta = 'θ ρ + (1 - θ) ρ0'
    ns.σ_ij = 'μ (utheta_i,j + utheta_j,i) - p δ_ij'  # diffusive term and pressure gradient

    ns.h = 1 / nelems
    ns.k = 'ρ h / μ'  # needs work, stabilization coeff

    # Residuals
    res = domain.integral('ρbasis_n ((ρ - ρ0) / dt + (ρtheta utheta_k)_,k) d:x' @ ns, degree=4)
    res += domain.integral('(ubasis_ni (((ρ u_i - ρ0 u0_i) / dt) + (ρtheta utheta_i utheta_j)_,j) + ubasis_ni,j σ_ij) d:x' @ ns, degree=4)
    # res += domain.integral('k ubasis_ni,j (u_j / sqrt(u_k u_k + 1e-12)) (ρ
    # (u_i - u0_i) / dt + (ρ u_i u_j - σ_ij)_,j) d:x' @ ns, degree=4) # SUPG
    # stabilization

    sqr = domain.boundary['right'].integral('(u_0 - 1)^2' @ ns, degree=4)
    cons0 = solver.optimize('lhs', sqr, droptol=1e-14)

    res0 = res
    lhs0 = np.zeros(res.shape)

    # Time loop
    t = 0.0
    timestep = 0
    bezier = domain.sample('bezier', 2)

    f = open("watchpoint.txt", "w")

    while participant.is_coupling_ongoing():

        p_read = participant.read_data(mesh_name, pressure_name, vertex_ids, precice_dt)
        treelog.info(f"p_read shape: {np.shape(p_read)}")
        treelog.info(f"p_read: {p_read}")
        ns.pin = p_read[0]
        res = res0 + domain.boundary['left'].integral('pin ubasis_ni n_i d:x' @ ns, degree=4)
        cons = min(0.2 * t, 1) * cons0

        if participant.requires_writing_checkpoint():
            lhscheckpoint = lhs0

        with treelog.context('stokes'):
            lhs = solver.newton('lhs', residual=res, constrain=cons, arguments=dict(lhs0=lhs0), lhs0=lhs0).solve(1e-08)
            treelog.info(f"lhs shape: {np.shape(lhs)}, type: {type(lhs)}")
            x, p, ρ, u = bezier.eval(['x_i', 'p', 'ρ', 'u_i'] @ ns, arguments=dict(lhs=lhs))

        treelog.info(f"velocity shape: {np.shape(u)}")
        treelog.info(f"velicity: {u}")
        write_vel = [[0, u[0][0], 0]]
        treelog.info(f"write_vel: {write_vel}")
        # write new velocities to 3D solver
        treelog.info(f"vertex ids: {vertex_ids}")
        participant.write_data(mesh_name, velocity_name, vertex_ids, write_vel)
        participant.advance(precice_dt)
        precice_dt = participant.get_max_time_step_size()

        if participant.requires_reading_checkpoint():
            lhs0 = lhscheckpoint
        else:
            lhs0 = lhs
            timestep += timestep
            x, p, ρ, u = bezier.eval(['x_i', 'p', 'ρ', 'u_i'] @ ns, lhs=lhs)
            f.write("%e; %e; %e; %e; %e; %e; %e\n" % (t, p[0], u[0], p[-1], u[-1], p[199], u[199]))
            f.flush()
            t += precice_dt
    participant.finalize()
    f.close()


if __name__ == '__main__':
    cli.run(main)
