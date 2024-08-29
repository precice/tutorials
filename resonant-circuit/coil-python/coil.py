import precice
import numpy as np
from scipy.integrate import solve_ivp

# Initialize and configure preCICE
participant = precice.Participant("Coil", "../precice-config.xml", 0, 1)

# Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
mesh_name = "Coil-Mesh"

dimensions = participant.get_mesh_dimensions(mesh_name)

vertex_ids = np.array([participant.set_mesh_vertex(mesh_name, np.zeros(dimensions))])

# Data IDs
read_data_name = "Voltage"
write_data_name = "Current"

# Simulation parameters and initial condition
C = 2                      # Capacitance of the capacitor
L = 1                      # Inductance of the coil
t0 = 0                     # Initial simulation time
Io = np.array([1])         # Initial current
phi = 0                    # Phase of the signal

w0 = 1 / np.sqrt(L * C)           # Resonant frequency
I0 = Io * np.cos(phi)           # Initial condition for I

# to estimate cost
global f_evals
f_evals = 0


def f_I(dt, max_allowed_dt):
    global f_evals
    f_evals += 1
    if dt > max_allowed_dt:  # read_data will fail, if dt is outside of window
        return np.nan

    U = participant.read_data(mesh_name, read_data_name, vertex_ids, dt)
    return U / L       # Time derivative of I; ODE determining capacitor


# Initialize simulation
if participant.requires_initial_data():
    participant.write_data(mesh_name, write_data_name, vertex_ids, I0)

participant.initialize()

solver_dt = participant.get_max_time_step_size()

# Start simulation
t = t0
while participant.is_coupling_ongoing():

    # Record checkpoint if necessary
    if participant.requires_writing_checkpoint():
        I0_checkpoint = I0
        t_checkpoint = t

    # Make Simulation Step
    precice_dt = participant.get_max_time_step_size()
    dt = min([precice_dt, solver_dt])
    t_span = [t, t + dt]
    sol = solve_ivp(lambda t, y: f_I(t - t_span[0], dt), t_span, I0, dense_output=True, rtol=1e-12, atol=1e-12)

    # Exchange data
    evals = max(len(sol.t), 3)  # at least do 3 substeps to allow cubic interpolation
    for i in range(evals):
        I0 = sol.sol(t_span[0] + (i + 1) * dt / evals)
        participant.write_data(mesh_name, write_data_name, vertex_ids, np.array(I0))
        participant.advance(dt / evals)

    t = t + dt

    # Recover checkpoint if not converged
    if participant.requires_reading_checkpoint():
        I0 = I0_checkpoint
        t = t_checkpoint

# Stop coupling
participant.finalize()


def I_analytical(t): return Io * np.cos(t * w0 + phi)


error = I0 - I_analytical(t)
print(f"{error=}")
print(f"{f_evals=}")
