import precice
import numpy as np
from scipy.integrate import solve_ivp

# Initialize and configure preCICE
participant = precice.Participant("ParticipantI", "../precice-config.xml", 0, 1)

# Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
mesh_name ="MeshI"

dimensions = participant.get_mesh_dimensions(mesh_name)

vertex_ids = np.array([participant.set_mesh_vertex(mesh_name, np.zeros(dimensions))])

# Data IDs
read_data_name = "U"
write_data_name = "I"

# Simulation parameters and initial condition
C = 2                      # Capacitance of the capacitor
L = 1                      # Inductance of the coil
t0 = 0                     # Initial simulation time
t_max = 10                 # End simulation time
Io = np.array([1])                     # Initial current
phi = 0                    # Phase of the signal

w0 = 1/np.sqrt(L*C)           # Resonant frequency
I0 = Io*np.cos(phi)           # Initial condition for I

def f_I(t_span, t):
    if t < t_span[1]:
        dt = t - t_span[0]
    else:
        dt = t_span[1] - t_span[0]

    U = participant.read_data(mesh_name, read_data_name, vertex_ids, dt)
    return U/L;       # Time derivative of I; ODE determining capacitor

# Initialize simulation
if participant.requires_initial_data():
    participant.write_data(mesh_name, write_data_name, vertex_ids, I0)

participant.initialize()

# Start simulation
t = t0
while participant.is_coupling_ongoing():

    # Record checkpoint if necessary
    if participant.requires_writing_checkpoint():
        I0_checkpoint = I0
        t_checkpoint = t

    # Make Simulation Step
    dt = participant.get_max_time_step_size()
    t_span = [t, t+dt]
    print(f"solving time {t_span}")
    sol = solve_ivp(lambda t, y: f_I(t_span, t), t_span, I0, dense_output=True,r_tol=1e-6, a_tol=1e-9)

    # Exchange data
    evals = 10
    for i in range(evals):
        I0 = sol.sol(t_span[0] + (i+1)*dt/evals)
        participant.write_data(mesh_name, write_data_name, vertex_ids, np.array(I0))
        participant.advance(dt/evals)

    t = t + dt

    # Recover checkpoint if not converged
    if participant.requires_reading_checkpoint():
        I0 = I0_checkpoint
        t = t_checkpoint

# Stop coupling
participant.finalize()

import math

assert(math.isclose(t, t_max))
I_analytical = lambda t: Io*np.cos(t*w0 + phi)
error = I0 - I_analytical(t_max)
print(f"error: {error}")