import numpy as np
import matplotlib.pyplot as plt
import treelog
from nutils import mesh, function, solver, cli
import csv

def main(nelems=200, dt=.005, tmax = 20.0, refdensity=1e3, refpressure=101325.0, psi=1e-6, viscosity=1e-3, theta=0.5):

    domain, geom = mesh.rectilinear([np.linspace(0,1000,nelems+1)])

    ns = function.Namespace()
    ns.x = geom
    ns.dt = dt
    ns.θ = theta
    ns.ρref = refdensity
    ns.pref = refpressure
    ns.pin = 98100
    ns.μ = viscosity
    ns.ψ = psi

    
    ns.ubasis, ns.ρbasis = function.chain([domain.basis('std', degree=2).vector(domain.ndims), domain.basis('std', degree=1)])

    ns.u_i = 'ubasis_ni ?lhs_n'
    ns.u0_i = 'ubasis_ni ?lhs0_n'
    ns.ρ = 'ρref + ρbasis_n ?lhs_n'
    ns.ρ0 = 'ρref + ρbasis_n ?lhs0_n'
    ns.p = 'pref + (ρ - ρref) / ψ'
    ns.utheta_i = 'θ u_i + (1 - θ) u0_i'
    ns.ρtheta = 'θ ρ + (1 - θ) ρ0'
    ns.σ_ij = 'μ (utheta_i,j + utheta_j,i) - p δ_ij' # diffusive term and pressure gradient

    ns.h = 1 / nelems
    ns.k = 'ρ h / μ' # needs work, stabilization coeff


    # Residuals
    res = domain.integral('ρbasis_n ((ρ - ρ0) / dt + (ρtheta utheta_k)_,k) d:x' @ ns, degree=4)
    res += domain.integral('(ubasis_ni (((ρ u_i - ρ0 u0_i) / dt) + (ρtheta utheta_i utheta_j)_,j) + ubasis_ni,j σ_ij) d:x' @ ns, degree=4)

    #res += domain.integral('ubasis_ni (ρtheta utheta_i utheta_j)_,j d:x' @ ns, degree=4)
    #res += domain.integral('ubasis_ni,j σ_ij d:x' @ ns, degree=4)
    res += domain.boundary['left'].integral('pin ubasis_ni n_i d:x' @ ns, degree=4)
    
    #res = domain.integral('ρbasis_n (δt(ρ) + t((ρ u_k)_,k)) d:x' @ ns, degree=4) # mass balance
    #res += domain.integral('(ubasis_ni (δt(ρ u_i) + t((ρ u_i u_j)_,j)) + ubasis_ni,j t(σ_ij)) d:x' @ ns, degree=4) # momentum balance
    #res += domain.boundary['left'].integral('pin ubasis_ni n_i d:x' @ ns, degree=4) # pressure set at inlet
    # res += domain.integral('k ubasis_ni,j (u_j / sqrt(u_k u_k + 1e-12)) (ρ (u_i - u0_i) / dt + (ρ u_i u_j - σ_ij)_,j) d:x' @ ns, degree=4) # SUPG stabilization


    sqr = domain.boundary['right'].integral('(u_0 - 1)^2' @ ns, degree=4)
    cons0 = solver.optimize('lhs', sqr, droptol=1e-14)

    lhs0 = np.zeros(res.shape)

    # Time loop
    t = 0.0
    timestep = 0
    next_plot_time = 0.0
    bezier = domain.sample('bezier', 2)
    results = []
    outlet_times = []          
    outlet_pressures = []      
    
    f = open("watchpoint.txt", "w")

    while t < tmax:
        treelog.info(f'Timestep {timestep}, Time = {t:.3f}')
        cons = min(0.2*t, 1)*cons0
        lhs = solver.newton('lhs', residual=res, constrain=cons, arguments=dict(lhs0=lhs0), lhs0=lhs0).solve(1e-08)
        lhs0 = lhs
        t += dt
        timestep += 1

        # Evaluate and store results
        x, p, ρ, u = bezier.eval(['x_i', 'p', 'ρ', 'u_i'] @ ns, lhs=lhs)
        f.write("%e; %e; %e\n" % (t, p[399], u[199]))
        f.flush()
        results.append((t, x[:, 0], u, p, ρ))
        outlet_times.append(t)                             
        outlet_pressures.append(p[-1])                     


        # Save plot at each 0.5s
        if t >= next_plot_time - 1e-8:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
            ax1.plot(x[:, 0], u, label='Velocity')
            ax2.plot(x[:, 0], p, label='Pressure', color='orange')
            ax1.set_ylabel('u')
            ax2.set_ylabel('p')
            ax2.set_xlabel('x')
            ax1.legend()
            ax2.legend()
            plt.tight_layout()
            plt.savefig(f'results/flow_t{t:.1f}.png')
            plt.close()
            next_plot_time += 0.5
    f.close()

    # Plot outlet pressure over the full simulation time
    plt.figure(figsize=(6, 4))
    plt.plot(outlet_times, outlet_pressures, label='Outlet Pressure')
    plt.xlabel('Time [s]')
    plt.ylabel('Outlet Pressure [Pa]')
    plt.title('Outlet Pressure vs Time (0–20s)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('results/outlet_pressure_vs_time_full.png')   # <-- NEW
    plt.close()

    treelog.info("Saved outlet pressure plot.")             # <-- NEW



if __name__ == '__main__':
    cli.run(main)