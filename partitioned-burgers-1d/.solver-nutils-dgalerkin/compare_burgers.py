import os
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))

data_dir_dgalerkin = os.path.join(DATA_DIR, "solver-nutils-dgalerkin", "data-training")
data_dir_diffuse = os.path.join(DATA_DIR, "solver-scipy-fvolumes", "data-training")

# file names are the same, folders are different
files_ = os.listdir(data_dir_dgalerkin)
files_ = sorted(f for f in files_ if f.endswith(".npz"))


file_num = 0
timestep = 0

file_dgalerkin = os.path.join(data_dir_dgalerkin, files_[file_num])
file_diffuse = os.path.join(data_dir_diffuse, files_[file_num])

data_dgalerkin = np.load(file_dgalerkin)["Solver-Mesh-1D-Internal"]
data_diffuse = np.load(file_diffuse)["Solver-Mesh-1D-Internal"]

print(f"DG data shape: {data_dgalerkin.shape}")
print(f"FV data shape: {data_diffuse.shape}")

plt.figure(figsize=(12, 5))
plt.plot(data_dgalerkin[timestep, :], label='DG Solver', color='blue')
plt.plot(data_diffuse[timestep, :], label='FV Solver', color='orange', linestyle='--')
plt.title(f'Comparison at Timestep {timestep}')
plt.xlabel('Spatial Position')
plt.ylabel('Solution Value')
plt.legend()
plt.savefig(os.path.join(DATA_DIR, f'comparison_timestep_{timestep}.png'))

# plot the imshow with unified colormap
vmin = min(data_dgalerkin.min(), data_diffuse.min())
vmax = max(data_dgalerkin.max(), data_diffuse.max())

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.imshow(data_dgalerkin.T, aspect='auto', cmap='viridis', vmin=vmin, vmax=vmax)
plt.title('DG Solver Evolution')
plt.ylabel('Spatial Position')
plt.xlabel('Timestep')
plt.colorbar(label='Solution Value')
plt.subplot(1, 2, 2)
plt.imshow(data_diffuse.T, aspect='auto', cmap='viridis', vmin=vmin, vmax=vmax)
plt.title('FV Solver Evolution')
plt.ylabel('Spatial Position')
plt.xlabel('Timestep')
plt.colorbar(label='Solution Value')
plt.tight_layout()
plt.show()
