import numpy as np
import json
import os
import argparse
import matplotlib.pyplot as plt

def _generate_initial_condition(x_coords, ic_config, epoch):
    np.random.seed(epoch)
    ic_values = np.zeros(len(x_coords))
    if ic_config["type"] == "sinusoidal":
        num_modes = ic_config.get("num_modes", 1)
        superpositions = np.random.randint(2, num_modes + 1)
        for _ in range(superpositions):
            amp = np.random.uniform(0.1, 2)
            k = np.random.randint(ic_config["wavenumber_range"][0], ic_config["wavenumber_range"][1] + 1)
            phase_shift = np.random.uniform(0, 2 * np.pi)
            ic_values += amp * np.sin(2 * np.pi * k * x_coords + phase_shift)
    return ic_values

def project_initial_condition(domain_min, domain_max, nelems, ic_config, epoch):
    # 1. Generate a high-resolution "truth" on a fine grid
    fine_res = nelems * 10
    fine_x = np.linspace(domain_min, domain_max, fine_res, endpoint=False)
    fine_u = _generate_initial_condition(fine_x, ic_config, epoch)

    # 2. Average the high-resolution truth over each coarse cell
    u_projected = np.zeros(nelems)
    for i in range(nelems):
        cell_start = i * 10
        cell_end = (i + 1) * 10
        u_projected[i] = np.mean(fine_u[cell_start:cell_end])
        
    return u_projected

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate initial conditions for the Burgers equation simulation.")
    parser.add_argument("--epoch", type=int, default=0, help="Seed for the random number generator to ensure reproducibility.")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(script_dir, "ic_params.json"), 'r') as f:
        config = json.load(f)
    
    ic_config = config["initial_conditions"]
    domain_config = config["domain"]

    # full domain
    full_domain_min = domain_config["full_domain_min"]
    full_domain_max = domain_config["full_domain_max"]
    nelems_total = domain_config["nelems_total"]

    # Generate IC
    initial_condition = project_initial_condition(full_domain_min, full_domain_max, nelems_total, ic_config, args.epoch)

    output_path = os.path.join(script_dir, "initial_condition.npz")
    np.savez(output_path, initial_condition=initial_condition)

    plt.figure(figsize=(8, 4))
    x_coords = np.linspace(full_domain_min, full_domain_max, nelems_total, endpoint=False)
    plt.figure(figsize=(10, 5))
    plt.plot(x_coords, initial_condition, marker='.', linestyle='-')
    plt.xlabel('Spatial Coordinate (x)')
    plt.ylabel('Solution Value (u)')
    plt.grid(True)
    plt.savefig(os.path.join(script_dir, "initial_condition.png"))
    print(f"Initial condition and plot saved to {output_path}")
    plt.close()
