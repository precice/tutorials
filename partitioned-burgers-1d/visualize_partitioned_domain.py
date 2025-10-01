import numpy as np
import matplotlib.pyplot as plt
import os

TIMESTEP_TO_PLOT = 50 #eg. 0, 1, ..., n, ... ,-1

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
DIRICHLET_DATA_PATH = os.path.join(CASE_DIR, "dirichlet-scipy", "dirichlet.npz")
NEUMANN_DATA_PATH = os.path.join(CASE_DIR, "neumann-scipy", "neumann.npz")
# NEUMANN_DATA_PATH = os.path.join(CASE_DIR, "surrogate-burgers", "surrogate.npz")

GROUND_TRUTH_DATA_PATH = os.path.join(CASE_DIR, "solver-scipy-fvolumes", "full_domain.npz")
if os.path.exists(GROUND_TRUTH_DATA_PATH):
    print(f"Found ground truth data at {GROUND_TRUTH_DATA_PATH}")
    gt_exists = True
else:
    print(f"Ground truth data not found at {GROUND_TRUTH_DATA_PATH}.\nPlease run python3 solver-scipy-fvolumes/solver.py None.")
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
    print(f"Loading ground truth data from {GROUND_TRUTH_DATA_PATH}")
    data_gt = np.load(GROUND_TRUTH_DATA_PATH)
    coords_gt = data_gt['internal_coordinates']
    solution_gt = data_gt['Solver-Mesh-1D-Internal']

    print(f"Ground truth shape: {solution_gt.shape}")

# --- plot single timestep ---
plt.figure(figsize=(10, 5), dpi=200)
plt.plot(full_coords, full_solution_history[TIMESTEP_TO_PLOT, :], marker='.', linestyle='-', label='Partitioned domain')
if gt_exists:
    plt.plot(coords_gt[:, 0], solution_gt[TIMESTEP_TO_PLOT, :], marker='x', linestyle='--', alpha=0.5, c="red", label='Ground truth')
plt.title(f'Solution at Timestep {TIMESTEP_TO_PLOT}')
plt.xlabel('Spatial Coordinate (x)')
plt.ylabel('Solution Value (u)')
plt.grid(True)
plt.legend()
plt.savefig(os.path.join(CASE_DIR, f'full_domain_timestep_slice.png'))
print(f"Saved plot to full_domain_timestep_slice.png")

if gt_exists:
    # residual
    residual = full_solution_history[TIMESTEP_TO_PLOT, :] - solution_gt[TIMESTEP_TO_PLOT, :]
    l2_norm_residual = np.linalg.norm(residual)
    l2_norm_gt = np.linalg.norm(solution_gt[TIMESTEP_TO_PLOT, :])
    relative_residual = l2_norm_residual / l2_norm_gt if l2_norm_gt > 1e-9 else 0.0

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
    print("Ground truth u at interface:")
    print(f"  t=0: u = {val_at_interface_t0:8.4f}, du/dx = {grad_at_interface_t0:8.4f}")
    print(f"  t={TIMESTEP_TO_PLOT}: u = {val_at_interface_plot:8.4f}, du/dx = {grad_at_interface_plot:8.4f}")
    print()
    print(f"Residual at t={TIMESTEP_TO_PLOT}:")
    print(f"  L2 Norm of Absolute Residual: {l2_norm_residual:10.6e}")
    print(f"  L2 Norm of Relative Residual: {relative_residual:10.6e}")
    print("---")

# --- plot gradient at single timestep ---
solution_slice = full_solution_history[TIMESTEP_TO_PLOT, :]
du_dx = np.gradient(solution_slice, full_coords)

plt.figure(figsize=(10, 5), dpi=200)
plt.plot(full_coords, du_dx, marker='.', linestyle='-', label='Partitioned domain')

if gt_exists:
    solution_gt_slice = solution_gt[TIMESTEP_TO_PLOT, :]
    du_dx_gt = np.gradient(solution_gt_slice, coords_gt[:, 0])
    plt.plot(coords_gt[:, 0], du_dx_gt, marker='x', linestyle='--', alpha=0.5, c="red", label='Ground truth')

plt.title(f'Gradient (du/dx) at Timestep {TIMESTEP_TO_PLOT}')
plt.xlabel('Spatial Coordinate (x)')
plt.ylabel('Gradient Value (du/dx)')
plt.grid(True)
plt.legend()
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
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, 'full_domain_evolution.png'))
print("Saved plot to full_domain_evolution.png")
plt.close()