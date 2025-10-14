import torch
import numpy as np
import precice
import os
import sys
import json

from model import CNN_RES
from config import INPUT_SIZE, OUTPUT_SIZE, HIDDEN_SIZE, NUM_RES_BLOCKS, KERNEL_SIZE, MODEL_NAME, ACTIVATION

GHOST_CELLS = INPUT_SIZE - OUTPUT_SIZE

def main():
    participant_name = "Neumann"
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    case_dir = os.path.abspath(os.path.join(script_dir, '..'))
    config_path = os.path.join(case_dir, "precice-config.xml")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = CNN_RES(
        hidden_channels=HIDDEN_SIZE,
        num_blocks=NUM_RES_BLOCKS,
        kernel_size=KERNEL_SIZE,
        activation=ACTIVATION,
        ghost_cells=GHOST_CELLS
    )

    if not os.path.exists(MODEL_NAME):
        print(f"Model file not found at {MODEL_NAME}")
        sys.exit(1)

    model.load_state_dict(torch.load(MODEL_NAME))
    model.to(device)
    model.eval()
    print("Neural surrogate model loaded successfully.")

    # Read initial condition
    with open(os.path.join(case_dir, "ic_params.json"), 'r') as f:
        domain_config = json.load(f)["domain"]

    nelems_total = domain_config["nelems_total"]
    nelems_local = nelems_total // 2
    full_domain_min = domain_config["full_domain_min"]
    full_domain_max = domain_config["full_domain_max"]
    dx = (full_domain_max - full_domain_min) / nelems_total

    ic_data = np.load(os.path.join(case_dir, "initial_condition.npz"))
    full_ic = ic_data['initial_condition']
    u = full_ic[nelems_local:]
    solution_history = {0.0: u.copy()}

    # Set domain and preCICE setup
    participant = precice.Participant(participant_name, config_path, 0, 1)

    mesh_name = "Neumann-Mesh"
    read_data_name = "Gradient"
    write_data_name = "Velocity"
    local_domain_min = full_domain_min + nelems_local * dx
    coupling_point = [[local_domain_min, 0.0]]
    vertex_id = participant.set_mesh_vertices(mesh_name, coupling_point)

    if participant.requires_initial_data():
        participant.write_data(mesh_name, write_data_name, vertex_id, [u[0]])

    participant.initialize()

    dt = participant.get_max_time_step_size()
    t = 0.0
    saved_t = 0.0
    t_index = 0
    solution_history = {int(0): u.copy()}

    # Main Coupling Loop 
    with torch.no_grad():
        while participant.is_coupling_ongoing():
            if participant.requires_writing_checkpoint():
                saved_u = u.copy()
                saved_t = t
                saved_t_index = t_index
            if participant.requires_reading_checkpoint():
                u = saved_u.copy()
                t = saved_t
                t_index = saved_t_index
                            
            du_dx_recv = participant.read_data(mesh_name, read_data_name, vertex_id, dt)[0]
            
            # Calculate ghost cell value from received gradient
            bc_left = u[0] - dx * du_dx_recv
            bc_right = u[-1]  # Zero gradient on right wall

            u_padded = np.empty(len(u) + 2)
            u_padded[0] = bc_left
            u_padded[-1] = bc_right
            u_padded[1:-1] = u
            
            input_tensor = torch.from_numpy(u_padded).float().unsqueeze(0).unsqueeze(0).to(device)
            
            output_tensor = model(input_tensor)
            u = output_tensor.squeeze().cpu().numpy()

            bc_left = u[0] - dx * du_dx_recv

            du_dx = (u[0] - bc_left) / dx
            u_interface = (bc_left + u[0]) / 2

            participant.write_data(mesh_name, write_data_name, vertex_id, [u[0]])

            print(f"[{participant_name:9s}] t={t:6.4f} | u_coupling={u_interface:8.4f} | du_dx={du_dx:8.4f}")

            t = saved_t + dt
            t_index += 1
            solution_history[t_index] = u.copy()
            participant.advance(dt)

    # Finalize and save data to npz array
    participant.finalize()

    run_dir = os.getcwd()
    output_filename = os.path.join(run_dir, "surrogate.npz")
    
    cell_centers_x = np.linspace(local_domain_min + dx/2, local_domain_min + (nelems_local - 0.5) * dx, nelems_local)
    internal_coords = np.array([cell_centers_x, np.zeros(nelems_local)]).T

    sorted_times_index = sorted(solution_history.keys())
    final_solution = np.array([solution_history[t_index] for t_index in sorted_times_index])

    np.savez(
        output_filename,
        internal_coordinates=internal_coords,
        **{"Solver-Mesh-1D-Internal": final_solution}
    )
    print(f"[Surrogate] Results saved to {output_filename}")

if __name__ == '__main__':
    main()