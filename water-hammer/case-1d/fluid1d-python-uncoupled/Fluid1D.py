import numpy as np
import matplotlib.pyplot as plt
import treelog
from nutils import mesh, function, solver, cli
import os


def main(nelems=200, dt=.005, tmax=20.0, refdensity=1e3, refpressure=101325.0, psi=1e-6, viscosity=1e-3, theta=0.55):
    """
    Simulates 1D unsteady compressible flow discretized using Nutils.
    Solves the momentum and mass equations with a theta time integration scheme.

    Parameters:
        nelems (int): Number of 1D elements in the domain.
        dt (float): Time step size.
        tmax (float): Final simulation time.
        refdensity (float): Reference fluid density.
        refpressure (float): Reference pressure.
        psi (float): Inverse of the speed of sound.
        viscosity (float): Dynamic viscosity of the fluid.
        theta (float): Time integration factor (theta-method; 0.5 = Crank-Nicolson).
    """

    # Create results directory
    os.makedirs('results', exist_ok=True)

    # Generate 1D mesh over domain [0, 1000] with `nelems` elements
    domain, geom = mesh.rectilinear([np.linspace(0, 1000, nelems + 1)])

    # Set up function namespace
    ns = function.Namespace()
    ns.x = geom
    ns.dt = dt
    ns.θ = theta
    ns.ρref = refdensity
    ns.pref = refpressure
    ns.pin = 98100             # Inlet pressure [Pa]
    ns.μ = viscosity
    ns.ψ = psi

    # Define basis functions: quadratic for velocity, linear for density
    ns.ubasis, ns.ρbasis = function.chain([
        domain.basis('std', degree=2).vector(domain.ndims),   # ubasis_ni
        domain.basis('std', degree=1)                         # ρbasis_n
    ])

    # Field definitions
    ns.u_i = 'ubasis_ni ?lhs_n'           # Current velocity field
    ns.u0_i = 'ubasis_ni ?lhs0_n'         # Previous velocity field
    ns.ρ = 'ρref + ρbasis_n ?lhs_n'       # Current density field
    ns.ρ0 = 'ρref + ρbasis_n ?lhs0_n'     # Previous density field
    ns.p = 'pref + (ρ - ρref) / ψ'        # Pressure via equation of state
    ns.utheta_i = 'θ u_i + (1 - θ) u0_i'  # Time-blended velocity (θ-scheme)
    ns.ρtheta = 'θ ρ + (1 - θ) ρ0'        # Time-blended density
    ns.σ_ij = 'μ (utheta_i,j + utheta_j,i) - p δ_ij'  # Stress tensor: viscous + pressure

    # === Residuals ===

    # Mass conservation residual
    res = domain.integral(
        'ρbasis_n ((ρ - ρ0) / dt + (ρtheta utheta_k)_,k) d:x' @ ns, degree=4)

    # Momentum conservation residual
    res += domain.integral(
        '(ubasis_ni (((ρ u_i - ρ0 u0_i) / dt) + (ρtheta utheta_i utheta_j)_,j) + ubasis_ni,j σ_ij) d:x' @ ns,
        degree=4)

    # Inlet boundary condition: impose pressure on left boundary
    res += domain.boundary['left'].integral(
        'pin ubasis_ni n_i d:x' @ ns, degree=4)

    # Outlet velocity constraint
    sqr = domain.boundary['right'].integral('(u_0 - 1)^2' @ ns, degree=4)
    cons0 = solver.optimize('lhs', sqr, droptol=1e-14)

    # Initial guess of lhs (left hand side)
    lhs0 = np.zeros(res.shape)

    # === Time loop setup ===
    t = 0.0
    timestep = 0
    next_plot_time = 0.0
    bezier = domain.sample('bezier', 2)  # Sample points for plotting

    # Open file to track watchpoint pressure and velocity at inlet/outlet/midpoint
    f = open("watchpoint.txt", "w")

    # === Time loop ===
    while t < tmax:
        treelog.info(f'Timestep {timestep}, Time = {t:.3f}')

        # Outlet velocity increases linearly from 0 to 1 m/s in 5 seconds, then remains 1 m/s
        cons = min(0.2 * t, 1) * cons0

        # Solve system with Newton
        with treelog.context('stokes'):
            lhs = solver.newton(
                'lhs', residual=res, constrain=cons,
                arguments=dict(lhs0=lhs0), lhs0=lhs0).solve(1e-08)

        # Update for next timestep
        lhs0 = lhs
        t += dt
        timestep += 1

        # Evaluate fields at bezier points
        x, p, ρ, u = bezier.eval(['x_i', 'p', 'ρ', 'u_i'] @ ns, lhs=lhs)

        # Log pressure and velocity at inlet, outlet, and midpoint
        # Format: time; p_in; u_in; p_out; u_out; p_mid; u_mid
        f.write("%e; %e; %e; %e; %e; %e; %e\n" % (t, p[0], u[0], p[-1], u[-1], p[199], u[199]))
        f.flush()

        # Save plot every 0.5s
        if t >= next_plot_time - 1e-8:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
            ax1.plot(x[:, 0], u, label='Velocity')
            ax2.plot(x[:, 0], p, label='Pressure', color='orange')
            ax1.set_ylabel('u [m/s]')
            ax2.set_ylabel('p [Pa]')
            ax2.set_xlabel('x [m]')
            ax1.set_title(f'Velocity at t = {t:.1f}s')
            ax2.set_title(f'Pressure at t = {t:.1f}s')
            ax1.legend()
            ax2.legend()
            plt.tight_layout()
            plt.savefig(f'results/flow_t{t:.1f}.png')
            plt.close()
            next_plot_time += 0.5

    # Close output file
    f.close()


if __name__ == '__main__':
    cli.run(main)
