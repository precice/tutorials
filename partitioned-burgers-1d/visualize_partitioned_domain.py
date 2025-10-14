import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('TIMESTEP_TO_PLOT', nargs='?', type=int, default=10, help="Timestep to plot, default is 10.")
parser.add_argument("--neumann", default="neumann-scipy/neumann.npz", help="Path to the neumann participant's data file relative to the case directory.")
args = parser.parse_args()
TIMESTEP_TO_PLOT = args.TIMESTEP_TO_PLOT

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
DIRICHLET_DATA_PATH = os.path.join(CASE_DIR, "dirichlet-scipy", "dirichlet.npz")
NEUMANN_DATA_PATH = os.path.join(CASE_DIR, args.neumann)

MONOLITHIC_DATA_PATH = os.path.join(CASE_DIR, "solver-scipy-fvolumes", "full_domain.npz")
if os.path.exists(MONOLITHIC_DATA_PATH):
    print(f"Found Monolithic data at {MONOLITHIC_DATA_PATH}")
    gt_exists = True
else:
    print(f"Monolithic data not found at {MONOLITHIC_DATA_PATH}.\nPlease run python3 solver-scipy-fvolumes/solver.py None.")
    gt_exists = False

print(f"Loading data from {DIRICHLET_DATA_PATH}")
data_d = np.load(DIRICHLET_DATA_PATH)
coords_d = data_d['internal_coordinates']
solution_d = data_d['Solver-Mesh-1D-Internal']

print(f"Loading data from {NEUMANN_DATA_PATH}")
data_n = np.load(NEUMANN_DATA_PATH)
coords_n = data_n['internal_coordinates']
solution_n = data_n['Solver-Mesh-1D-Internal']

full_coords = np.concatenate((coords_d[:, 0], coords_n[:, 0]))
full_solution_history = np.concatenate((solution_d, solution_n), axis=1)

print(f"Full domain shape: {full_solution_history.shape}")

if gt_exists:
    print(f"Loading Monolithic data from {MONOLITHIC_DATA_PATH}")
    data_gt = np.load(MONOLITHIC_DATA_PATH)
    coords_gt = data_gt['internal_coordinates']
    solution_gt = data_gt['Solver-Mesh-1D-Internal']

    print(f"Monolithic shape: {solution_gt.shape}")

# --- plot single timestep ---
plt.figure(figsize=(10, 5), dpi=300)
plt.plot(full_coords, full_solution_history[TIMESTEP_TO_PLOT, :], marker='.', markersize=8, linestyle='-', label=f'Partitioned Solution\n{args.neumann.split("/")[0]}')

if gt_exists:
    plt.plot(coords_gt[:, 0], solution_gt[TIMESTEP_TO_PLOT, :], marker='+', linestyle=':',  c="crimson", alpha=1, label='Monolithic Solution')
plt.title(f'Solution at Timestep {TIMESTEP_TO_PLOT}')
plt.xlabel('Spatial Coordinate (x)')
plt.ylabel('Solution Value (u)')
plt.grid(True)

u_max = np.max(full_solution_history[TIMESTEP_TO_PLOT, :])
u_min = np.min(full_solution_history[TIMESTEP_TO_PLOT, :])
u_interface = full_solution_history[TIMESTEP_TO_PLOT, full_coords.searchsorted(1.0)]
u_offset = (u_max-u_min)*0.125
plt.vlines(x=1, ymin=u_min-u_offset, ymax=u_interface-u_offset, color='gray', linestyle='--', label='Interface')
plt.vlines(x=1, ymin=u_interface+u_offset, ymax=u_max+u_offset*2, color='gray', linestyle='--')

plt.legend(loc='upper left')
plt.savefig(os.path.join(CASE_DIR, f'full_domain_timestep_slice.png'))
print(f"Saved plot to full_domain_timestep_slice.png")

if gt_exists:
    # residual
    residual = full_solution_history[TIMESTEP_TO_PLOT, :] - solution_gt[TIMESTEP_TO_PLOT, :]
    mse = np.mean(np.square(residual))
    mse_gt_vs_zero = np.mean(np.square(solution_gt[TIMESTEP_TO_PLOT, :]))
    relative_mse = mse / mse_gt_vs_zero if mse_gt_vs_zero > 1e-9 else 0.0

    nelems_total = solution_gt.shape[1]
    interface_idx = nelems_total // 2 - 1
    dx = coords_gt[1, 0] - coords_gt[0, 0]

    # t = 0
    u0_gt = solution_gt[0, :]
    val_at_interface_t0 = (u0_gt[interface_idx] + u0_gt[interface_idx + 1]) / 2.0
    grad_at_interface_t0 = (u0_gt[interface_idx + 1] - u0_gt[interface_idx]) / dx

    # t = TIMESTEP_TO_PLOT
    u_plot_gt = solution_gt[TIMESTEP_TO_PLOT, :]
    val_at_interface_plot = (u_plot_gt[interface_idx] + u_plot_gt[interface_idx + 1]) / 2.0
    grad_at_interface_plot = (u_plot_gt[interface_idx + 1] - u_plot_gt[interface_idx]) / dx

    print("---")
    print("Monolithic u at interface:")
    print(f"  t=0: u = {val_at_interface_t0:8.4f}, du/dx = {grad_at_interface_t0:8.4f}")
    print(f"  t={TIMESTEP_TO_PLOT}: u = {val_at_interface_plot:8.4f}, du/dx = {grad_at_interface_plot:8.4f}")
    print()
    print(f"Residual at t={TIMESTEP_TO_PLOT}:")
    print(f"  Mean Squared Error (MSE): {mse:10.6e}")
    print(f"  Relative MSE: {relative_mse:10.6e}")
    print("---")

# --- plot gradient at single timestep ---
solution_slice = full_solution_history[TIMESTEP_TO_PLOT, :]
du_dx = np.gradient(solution_slice, full_coords)

plt.figure(figsize=(10, 5), dpi=300)
plt.plot(full_coords, du_dx, marker='.',  markersize=8, linestyle='-', label=f'Partitioned Solution\n{args.neumann.split("/")[0]}')

if gt_exists:
    solution_gt_slice = solution_gt[TIMESTEP_TO_PLOT, :]
    du_dx_gt = np.gradient(solution_gt_slice, coords_gt[:, 0])
    plt.plot(coords_gt[:, 0], du_dx_gt, marker='+', linestyle=':',  c="crimson", alpha=1, label='Monolithic Solution')

plt.title(f'Gradient (du/dx) at Timestep {TIMESTEP_TO_PLOT}')
plt.xlabel('Spatial Coordinate (x)')
plt.ylabel('Gradient Value (du/dx)')
plt.grid(True)

u_max = np.max(du_dx)
u_min = np.min(du_dx)
u_interface = du_dx[full_coords.searchsorted(1.0)]
u_offset = (u_max-u_min)*0.125
plt.vlines(x=1, ymin=u_min-u_offset, ymax=u_interface-u_offset, color='gray', linestyle='--', label='Interface')
plt.vlines(x=1, ymin=u_interface+u_offset, ymax=u_max+u_offset, color='gray', linestyle='--')

plt.legend(loc='upper left')
plt.savefig(os.path.join(CASE_DIR, f'gradient_timestep_slice.png'))
print(f"Saved plot to gradient_timestep_slice.png")
plt.close()

# --- plot time evolution ---
plt.figure(figsize=(10, 6))
plt.imshow(full_solution_history.T, aspect='auto', cmap='viridis', origin='lower',
           extent=[0, full_solution_history.shape[0], full_coords.min(), full_coords.max()])
plt.colorbar(label='Solution Value (u)')
plt.title('Time Evolution of Partitioned Burgers eq.')
plt.xlabel('Timestep')
plt.ylabel('Spatial Coordinate (x)')
plt.xticks(np.arange(0, full_solution_history.shape[0], step=max(1, full_solution_history.shape[0]//10)))
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, 'full_domain_evolution.png'))
print("Saved plot to full_domain_evolution.png")
plt.close()