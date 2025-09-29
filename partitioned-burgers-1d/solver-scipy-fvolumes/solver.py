"""
This script solves the 1D viscous Burgers' equation for a partitioned domain.
- 'Dirichlet': Solves the left half of the domain.
- 'Neumann':   Solves the right half of the domain.

# |<-----------------Dirichlet-------------------->|<--------------------Neumann-------------------->|
# | bc_left = 0|  |  ...  |   | u[-1] <|> bc_right | bc_left <|> u[0] |   |  ...  |   | bc_right = 0 |
# |<------------------ u ------------------------->|<-------------------- u ------------------------>|
# |                              <u_interface>     |    <u_interface>                                |
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
import precice
import json
import os
import argparse
import time

def lax_friedrichs_flux(u_left, u_right): # Flux at cell interface between i-1/2 and i+1/2
    return 0.5 * (0.5 * u_left**2 + 0.5 * u_right**2) - 0.5 * (u_right - u_left)

def flux_function(u_left, u_right):
    return lax_friedrichs_flux(u_left, u_right)

def burgers_rhs(t, u, dx, C, bc_left, bc_right):
    # Ghost cells for BCs
    u_padded = np.empty(len(u) + 2)
    u_padded[0] = bc_left
    u_padded[-1] = bc_right
    u_padded[1:-1] = u

    flux = np.empty(len(u) + 1)
    for i in range(len(flux)):
        flux[i] = flux_function(u_padded[i], u_padded[i+1])

    # viscosity
    viscosity = C * (u_padded[2:] - 2 * u_padded[1:-1] + u_padded[:-2]) / dx**2
    return -(flux[1:] - flux[:-1]) / dx + viscosity

# def burgers_jacobian(t, u, dx, C, bc_left, bc_right):
#     n = len(u)
#     # Derivatives of the Lax-Friedrichs flux
#     # dF/du_L = 0.5 * u_L + 0.5
#     # dF/du_R = 0.5 * u_R - 0.5
#     
#     # Contribution of flux_{i+1/2} to the diagonal at i
#     d_flux_di_plus = -(0.5 * u + 0.5) / dx
#     # Contribution of flux_{i-1/2} to the diagonal at i
#     d_flux_di_minus = (0.5 * u - 0.5) / dx
#     
#     # Off-diagonals
#     d_flux_d_up1 = -(0.5 * u[1:] - 0.5) / dx
#     d_flux_d_um1 = (0.5 * u[:-1] + 0.5) / dx
#
#     # --- Viscous Part ---
#     d_visc_di = -2 * C / dx**2
#     d_visc_off = C / dx**2
#
#     main_diag = d_flux_di_plus + d_flux_di_minus + d_visc_di
#     upper_diag = d_flux_d_up1 + d_visc_off
#     lower_diag = d_flux_d_um1 + d_visc_off
#
#     # Construct Jacobian
#     jac = diags([main_diag, upper_diag, lower_diag], [0, 1, -1], shape=(n, n), format='csc')
#
#     return jac

def main(participant_name: str):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.abspath(os.path.join(script_dir, '..'))
    run_dir = os.getcwd()


    if participant_name == 'None':
        print("Participant not specified. Running full domain without preCICE")
        participant_name = None
    else:
        config_path = os.path.join(case_dir, "precice-config.xml")
        participant = precice.Participant(participant_name, config_path, 0, 1)

    # Read initial condition
    with open(os.path.join(case_dir, "ic_params.json"), 'r') as f:
        domain_config = json.load(f)["domain"]

    nelems_total = domain_config["nelems_total"]
    nelems_local = nelems_total // 2
    full_domain_min = domain_config["full_domain_min"]
    full_domain_max = domain_config["full_domain_max"]
    dx = (full_domain_max - full_domain_min) / nelems_total

    # Set domain and preCICE setup
    if participant_name == "Dirichlet":
        mesh_name = "Dirichlet-Mesh"
        read_data_name = "Velocity"
        write_data_name = "Gradient"
        local_domain_min = full_domain_min
        local_domain_max = full_domain_min + nelems_local * dx
        coupling_point = [[local_domain_max, 0.0]]
    elif participant_name == "Neumann":
        mesh_name = "Neumann-Mesh"
        read_data_name = "Gradient"
        write_data_name = "Velocity"
        local_domain_min = full_domain_min + nelems_local * dx
        local_domain_max = full_domain_max
        coupling_point = [[local_domain_min, 0.0]]
    else: #full domain run
        local_domain_min = full_domain_min
        local_domain_max = full_domain_max
        nelems_local = nelems_total
        dt = 0.01 # Fixed time step for standalone run
        t_end = 0.1

    if participant_name is not None:
        vertex_id = participant.set_mesh_vertices(mesh_name, coupling_point)
        participant.initialize()
        dt = participant.get_max_time_step_size()

    ic_data = np.load(os.path.join(case_dir, "initial_condition.npz"))
    full_ic = ic_data['initial_condition']
    if participant_name == "Dirichlet":
        u = full_ic[:nelems_local]
    elif participant_name == "Neumann":
        u = full_ic[nelems_local:]
    else:
        u = full_ic

    solution_history = {int(0): u.copy()}

    t = 0.0
    t_index = 0
    saved_t = 0.0
    C_viscosity = 1e-12
    aborted = False

    # --- Serial Coupling Loop ---
    if participant_name == "Dirichlet":
        while participant.is_coupling_ongoing():
            if participant.requires_writing_checkpoint():
                # print(f"[Dirichlet] Writing checkpoint at t={t:.4f}")
                saved_u = u.copy()
                saved_t = t
            if participant.requires_reading_checkpoint():
                u = saved_u.copy()
                t = saved_t
                # print(f"[Dirichlet] Reading checkpoint at t={t:.4f}")

            u_from_neumann = participant.read_data(mesh_name, read_data_name, vertex_id, dt)[0]

            bc_right = u_from_neumann

            t_end = t + dt
            solver_args = (dx, C_viscosity, 0, bc_right)
            sol = solve_ivp(burgers_rhs, (t, t_end), u, args=solver_args, method='BDF', t_eval=[t_end])
            u = sol.y[:, -1]
            # u = u + dt * burgers_rhs(t, u, *solver_args)


            gradient_to_send = (u_from_neumann - u[-1]) / dx
            flux_across_interface = flux_function(u[-1], bc_right)
     
            participant.write_data(mesh_name, write_data_name, vertex_id, [gradient_to_send])

            print(f"[{participant_name:9s}] t={t:6.4f} | u_coupling={u_from_neumann:8.4f} | grad_sent={gradient_to_send:8.4f} | flux_across={flux_across_interface:8.4f}")

            t = saved_t + dt
            t_index = int(t/dt)
            solution_history[t_index] = u.copy()
            participant.advance(dt)

    elif participant_name == "Neumann":
        while participant.is_coupling_ongoing():
            if participant.requires_writing_checkpoint():
                # print(f"[Neumann] Writing checkpoint at t={t:.4f}")
                saved_u = u.copy()
                saved_t = t
            if participant.requires_reading_checkpoint():
                u = saved_u.copy()
                t = saved_t
                # print(f"[Neumann] Reading checkpoint at t={t:.4f}")
                            
            du_dx_recv = participant.read_data(mesh_name, read_data_name, vertex_id, dt)[0]
            
            bc_left = u[0] - du_dx_recv * dx
            
            t_end = t + dt
            solver_args = (dx, C_viscosity, bc_left, 0)
            sol = solve_ivp(burgers_rhs, (t, t_end), u, args=solver_args, method='BDF', t_eval=[t_end])
            u = sol.y[:, -1]
            # u = u + dt * burgers_rhs(t, u, *solver_args)


            bc_left = u[0] - du_dx_recv * dx # Update bc_left to be consistent with the new state
            flux_across_interface = flux_function(bc_left, u[0])

            # if u[0] < 0 and t > 0:
            #     print(f"[{participant_name}] Only upwind scheme implemented, aborting!")
            #     aborted = True
            #     break

            participant.write_data(mesh_name, write_data_name, vertex_id, [u[0]])

            print(f"[{participant_name:9s}] t={t:6.4f} | u_coupling={u[0]:8.4f} | grad_recv={du_dx_recv:8.4f} | flux_across={flux_across_interface:8.4f}")

            t = saved_t + dt
            t_index = int(t/dt)
            solution_history[t_index] = u.copy()
            participant.advance(dt)

    if participant_name is not None:
        # Finalize and save data to npz array
        participant.finalize()
        output_filename = os.path.join(run_dir, f"{participant_name.lower()}.npz")
    else:
        output_filename = os.path.join(script_dir, "full_domain.npz")
        print("Starting standalone simulation without preCICE")
        bc_left, bc_right = 0, 0

        while t + dt < t_end:
            step_end = t + dt
            solver_args = (dx, C_viscosity, bc_left, bc_right)
            sol = solve_ivp(burgers_rhs, (t, step_end), u, args=solver_args, method='BDF', t_eval=[step_end], )
            u = sol.y[:, -1]
            # u = u + dt * burgers_rhs(t, u, *solver_args)
            
            t = t + dt
            t_index = int(t/dt)
            solution_history[t_index] = u.copy()
            print(f"[Standalone ] t={t:6.4f}")

    if not aborted:

        cell_centers_x = np.linspace(local_domain_min + dx/2, local_domain_max - dx/2, nelems_local)
        internal_coords = np.array([cell_centers_x, np.zeros(nelems_local)]).T

        sorted_times_index = sorted(solution_history.keys())
        final_solution = np.array([solution_history[t_index] for t_index in sorted_times_index])

        np.savez(
            output_filename,
            internal_coordinates=internal_coords,
            **{"Solver-Mesh-1D-Internal": final_solution}
        )
        print(f"Results saved to {output_filename}")
    else:
        raise RuntimeError("Simulation aborted.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("participant", help="Name of the participant", choices=['Dirichlet', 'Neumann', 'None'])
    args = parser.parse_args()
    main(args.participant)