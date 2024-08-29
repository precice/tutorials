import precice
import numpy as np
from scipy.integrate import solve_ivp

# Initialize and configure preCICE
participant = precice.Participant("Capacitor", "../precice-config.xml", 0, 1)

# Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
mesh_name = "Capacitor-Mesh"

dimensions = participant.get_mesh_dimensions(mesh_name)

vertex_ids = np.array([participant.set_mesh_vertex(mesh_name, np.zeros(dimensions))])

# Data IDs
read_data_name = "Current"
write_data_name = "Voltage"

# Simulation parameters and initial condition
C = 2                      # Capacitance of the capacitor
L = 1                      # Inductance of the coil
t0 = 0                     # Initial simulation time
t_max = 10                 # End simulation time
Io = np.array([1])         # Initial current
phi = 0                    # Phase of the signal

w0 = 1 / np.sqrt(L * C)           # Resonant frequency
U0 = -w0 * L * Io * np.sin(phi)     # Initial condition for U


def f_U(dt, max_allowed_dt):
    if dt > max_allowed_dt:  # read_data will fail, if dt is outside of window
        return np.nan

    I = participant.read_data(mesh_name, read_data_name, vertex_ids, dt)
    return -I / C      # Time derivative of U


# Initialize simulation
if participant.requires_initial_data():
    participant.write_data(mesh_name, write_data_name, vertex_ids, U0)

participant.initialize()

solver_dt = participant.get_max_time_step_size()

# Start simulation
t = t0
while participant.is_coupling_ongoing():

    # Record checkpoint if necessary
    if participant.requires_writing_checkpoint():
        U0_checkpoint = U0
        t_checkpoint = t

    # Make Simulation Step
    precice_dt = participant.get_max_time_step_size()
    dt = min([precice_dt, solver_dt])
    t_span = [t, t + dt]
    sol = solve_ivp(lambda t, y: f_U(t - t_span[0], dt), t_span, U0, dense_output=True, rtol=1e-12, atol=1e-12)

    # Exchange data
    evals = max(len(sol.t), 3)  # at least do 3 substeps to allow cubic interpolation
    for i in range(evals):
        U0 = sol.sol(t_span[0] + (i + 1) * dt / evals)
        participant.write_data(mesh_name, write_data_name, vertex_ids, np.array(U0))
        participant.advance(dt / evals)

    t = t + dt

    # Recover checkpoint if not converged
    if participant.requires_reading_checkpoint():
        U0 = U0_checkpoint
        t = t_checkpoint

# Stop coupling
participant.finalize()
