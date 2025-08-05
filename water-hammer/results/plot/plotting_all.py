import numpy as np
import matplotlib.pyplot as plt

def load_and_filter_1d(filename, tmin=6.8, tmax=8.2):
    data = np.loadtxt(filename, delimiter=';')
    time = data[:, 0]
    pressure = data[:, 1]
    mask = (time >= tmin) & (time <= tmax)
    return time[mask], pressure[mask]

def load_and_filter_1d_1d(filename, tmin=6.8, tmax=8.2):
    data = np.loadtxt(filename, delimiter=';')
    time = data[:, 0]
    pressure = data[:, 3]
    mask = (time >= tmin) & (time <= tmax)
    return time[mask], pressure[mask]

def load_and_filter_openfoam(filepath, tmin=6.8, tmax=8.2):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # Skip comment lines
    data = [line.strip().split() for line in lines if not line.startswith('#')]
    data = np.array(data, dtype=float)
    time = data[:, 0]
    value = data[:, 1]
    mask = (time >= tmin) & (time <= tmax)
    return time[mask], value[mask]

# Load and filter each file
t1, p1 = load_and_filter_openfoam('../I/p')                   # 1D-3D
t2, p2 = load_and_filter_1d('../II/watchpoint.txt')           # 1D, θ = 0.5
t3, p3 = load_and_filter_1d('../III/watchpoint.txt')          # 1D, θ = 0.55
t4, p4 = load_and_filter_1d('../IV/watchpoint.txt')           # 1D, θ = 1  
t5, p5 = load_and_filter_openfoam('../V/p')                   # 3D
t6, p6 = load_and_filter_1d_1d('../VI/watchpointright.txt')   # 1D-1D
t7, p7 = load_and_filter_openfoam('../VII/pright')            # 3D-3D

# Plot
plt.figure(figsize=(9, 5))
plt.plot(t1, p1, label='1D-3D Coupling', color='purple', linewidth=0.8)
plt.plot(t2, p2, label='1D; θ = 0.5', color='blue', linewidth=0.8)
plt.plot(t3, p3, label='1D; θ = 0.55', color='green', linewidth=0.8)
plt.plot(t4, p4, label='1D; θ = 1.0', color='red', linewidth=0.8)
plt.plot(t5, p5, label='3D', color='orange', linewidth=0.8)
plt.plot(t6, p6, label='1D-1D Coupling', color='brown', linewidth=0.8)
plt.plot(t7, p7, label='3D-3D Coupling', color='black', linewidth=0.8)

plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
plt.title('Outlet Pressure (6.8s - 8.2s)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('pout_timewindow.png')
plt.show()