"""
This script solves the 1D viscous Burgers' equation for a partitioned domain problem.
The two participants:
- 'Dirichlet': Solves the left half of the domain. Receives Dirichlet BC from 'Neumann'. Provides Neumann BC to 'Neumann'.
- 'Neumann':   Solves the right half of the domain.

# |<---------------------Dirichlet---------------------->|<--------------------Neumann----------------------->|
# | du_dx=0 |  ...  |  ...  |  ...  | u[-1] <|> bc_right | bc_left <|> u[0] |  ...  |  ...  |  ...  | du_dx=0 |
# |<----------------------- u -------------------------->|<----------------------- u ------------------------>|
# | bc_left |                          <u_interface>     |    <u_interface>                        | bc_right |
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, identity
from scipy.optimize import root
import precice
import json
import os
import argparse

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

def burgers_jacobian(t, u, dx, C, bc_left, bc_right):
    n = len(u)
    
    u_padded = np.empty(n + 2)
    u_padded[0] = bc_left
    u_padded[-1] = bc_right
    u_padded[1:-1] = u
    
    # --- viscosity (constant) ---
    d_visc_di = -2 * C / dx**2
    d_visc_off = C / dx**2

    # --- flux (Lax-Friedrichs) ---
    df_du_di = -1 / dx

    # Upper diagonal
    df_du_upper = -(0.5 * u_padded[2:n+1] - 0.5) / dx

    # Lower diagonal
    df_du_lower = (0.5 * u_padded[0:n-1] + 0.5) / dx

    main_diag = df_du_di + d_visc_di
    upper_diag = df_du_upper + d_visc_off
    lower_diag = df_du_lower + d_visc_off

    jac = diags([main_diag, upper_diag, lower_diag], [0, 1, -1], shape=(n, n), format='csc')

    return jac

def burgers_residual(u_new, u_old, dt, dx, C, bc_left, bc_right):
    return u_new - u_old - dt * burgers_rhs(0, u_new, dx, C, bc_left, bc_right)

def burgers_jacobian_residual(u_new, u_old, dt, dx, C, bc_left, bc_right):
    n = len(u_new)
    I = identity(n, format='csc')
    J_rhs = burgers_jacobian(0, u_new, dx, C, bc_left, bc_right)
    return (I - dt * J_rhs).toarray() # doesn't work with sparse matrix in root solver for some reason

class BoundaryWrapper:
    """
    Wrap the RHS and Jacobian to dynamically set BCs during the solve iterations with the updated state u. 
    """
    def __init__(self, dx, C, participant_name, u_from_neumann=None, du_dx_recv=None):
        self.dx = dx
        self.C = C
        self.participant_name = participant_name
        self.u_from_neumann = u_from_neumann
        self.du_dx_recv = du_dx_recv

    def bc_left(self, u):
        if self.participant_name == "Neumann":
            return u[0] - self.du_dx_recv * self.dx
        # zero gradient at outer boundary
        elif self.participant_name == "Dirichlet":
            return u[0]
        else:
            return u[0]
    
    def bc_right(self, u):
        if self.participant_name == "Dirichlet":
            return self.u_from_neumann
        # zero gradient at outer boundary
        elif self.participant_name == "Neumann":
            return u[-1]
        else:
            return u[-1]

    def rhs(self, t, u):
        bc_left = self.bc_left(u)
        bc_right = self.bc_right(u)
        return burgers_rhs(t, u, self.dx, self.C, bc_left, bc_right)

    def jac(self, t, u):
        bc_left = self.bc_left(u)
        bc_right = self.bc_right(u)
        
        J_rhs = burgers_jacobian(t, u, self.dx, self.C, bc_left, bc_right)
        return J_rhs
    
    def rhs_residual(self, u_new, u_old, dt):
        bc_left = self.bc_left(u_new)
        bc_right = self.bc_right(u_new)
        return burgers_residual(u_new, u_old, dt, self.dx, self.C, bc_left, bc_right)
    
    def jac_residual(self, u_new, u_old, dt):
        bc_left = self.bc_left(u_new)
        bc_right = self.bc_right(u_new)
        return burgers_jacobian_residual(u_new, u_old, dt, self.dx, self.C, bc_left, bc_right)

def main(participant_name: str):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.abspath(os.path.join(script_dir, '..'))
    run_dir = os.getcwd()

    config_path = os.path.join(case_dir, "precice-config.xml")

    if participant_name == 'None':
        # read precice config to get t_final and dt
        import re
        print("Participant not specified. Running full domain without preCICE")
        participant_name = None

        with open(config_path, 'r') as f:
            precice_config = f.read()
        max_time_match = re.search(r'<max-time\s+value="([^"]+)"\s*/>', precice_config)
        t_final = float(max_time_match.group(1))
        dt_match = re.search(r'<time-window-size\s+value="([^"]+)"\s*/>', precice_config)
        dt = float(dt_match.group(1))
        print(f"t_final = {t_final}, dt = {dt}")
    else:
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

    ic_data = np.load(os.path.join(case_dir, "initial_condition.npz"))
    full_ic = ic_data['initial_condition']
    if participant_name == "Dirichlet":
        u = full_ic[:nelems_local]
    elif participant_name == "Neumann":
        u = full_ic[nelems_local:]
    else:
        u = full_ic

    if participant_name is not None:
        vertex_id = participant.set_mesh_vertices(mesh_name, coupling_point)

        if participant.requires_initial_data():
            if participant_name == "Dirichlet":
                du_dx_send = (u[-1] - u[-2]) / dx # take forward difference inside domain for initial send
                participant.write_data(mesh_name, write_data_name, vertex_id, [du_dx_send])
            if participant_name == "Neumann":
                participant.write_data(mesh_name, write_data_name, vertex_id, [u[0]])

        participant.initialize()
        dt = participant.get_max_time_step_size()

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

            t_end = t + dt
            wrapper = BoundaryWrapper(dx, C_viscosity, "Dirichlet", u_from_neumann=u_from_neumann)

            # sol = solve_ivp(wrapper.rhs, (t, t_end), u, method='BDF', t_eval=[t_end], jac=wrapper.jac) # BDF higher order implicit timestepping
            # u = sol.y[:, -1]
            # u = u + dt * burgers_rhs(t, u, dx,  C_viscosity, wrapper.bc_left(u), wrapper.bc_right(u)) # Explicit Euler

            sol = root(wrapper.rhs_residual, u, args=(u, dt), jac=wrapper.jac_residual, method='hybr')
            u = sol.x

            bc_right = wrapper.bc_right(u)

            du_dx_send = (bc_right - u[-1]) / dx
            flux_across_interface = flux_function(u[-1], bc_right)
            u_interface = (u[-1] + bc_right) / 2
     
            participant.write_data(mesh_name, write_data_name, vertex_id, [du_dx_send])

            print(f"[{participant_name:9s}] t={t:6.4f} | u_coupling={u_interface:8.4f} | du_dx={du_dx_send:8.4f} | flux_across={flux_across_interface:8.4f}")

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
            
            t_end = t + dt
            wrapper = BoundaryWrapper(dx, C_viscosity, "Neumann", du_dx_recv=du_dx_recv)
            # sol = solve_ivp(wrapper.rhs, (t, t_end), u, method='BDF', t_eval=[t_end], jac=wrapper.jac) # BDF higher order implicit timestepping
            # u = sol.y[:, -1]
            # u = u + dt * burgers_rhs(t, u, dx,  C_viscosity, wrapper.bc_left(u), wrapper.bc_right(u)) # Explicit Euler

            sol = root(wrapper.rhs_residual, u, args=(u, dt), jac=wrapper.jac_residual, method='hybr')
            u = sol.x

            bc_left = wrapper.bc_left(u)
            flux_across_interface = flux_function(bc_left, u[0])
            du_dx = (u[0] - bc_left) / dx
            u_interface = (bc_left + u[0]) / 2

            participant.write_data(mesh_name, write_data_name, vertex_id, [u[0]])

            print(f"[{participant_name:9s}] t={t:6.4f} | u_coupling={u_interface:8.4f} | du_dx={du_dx:8.4f} | flux_across={flux_across_interface:8.4f}")

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

        while t + dt < t_final:

            print(f"[Standalone ] t={t:6.4f}")
            t_end = t + dt
            wrapper = BoundaryWrapper(dx, C_viscosity, "None")
            # sol = solve_ivp(wrapper.rhs, (t, t_end), u, method='BDF', t_eval=[t_end], jac=wrapper.jac) # BDF higher order implicit timestepping
            # u = sol.y[:, -1]
            # u = u + dt * burgers_rhs(t, u, dx,  C_viscosity, wrapper.bc_left(u), wrapper.bc_right(u)) # Explicit Euler

            sol = root(wrapper.rhs_residual, u, args=(u, dt), jac=wrapper.jac_residual, method='hybr')
            u = sol.x
            
            t = t + dt
            t_index = int(t/dt)
            solution_history[t_index] = u.copy()

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