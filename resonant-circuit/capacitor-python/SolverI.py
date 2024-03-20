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
Io = 1                     # Initial current
phi = 0                    # Phase of the signal

w0 = 1/np.sqrt(L*C)           # Resonant frequency
I0 = Io*np.cos(phi)           # Initial condition for I
U0 = -w0*L*Io*np.sin(phi)     # Initial condition for U

f_I = lambda t, I, U: U/L;       # Time derivative of I; ODE determining capacitor

# Initialize simulation
if participant.requires_initial_data():
    participant.write_data(mesh_name, write_data_name, vertex_ids, np.array([I0]))

participant.initialize()

dt = participant.get_max_time_step_size()

# Start simulation
t = t0 + dt
while participant.is_coupling_ongoing():

    # Record checkpoint if necessary
    if participant.requires_writing_checkpoint():
        I0_checkpoint = I0
        U0_checkpoint = U0

    # Make Simulation Step
    sol = solve_ivp(lambda t, y: f_I(t, y, U0), [t0, t], np.array([I0]))
    I0 = sol.y[-1,-1]
    print(sol.y)
    print(I0)

    # Exchange data
    participant.write_data(mesh_name, write_data_name, vertex_ids, np.array([I0]))
    participant.advance(dt)

    # Recover checkpoint if not converged, else finish time step
    if participant.requires_reading_checkpoint():
        I0 = I0_checkpoint
        U0 = U0_checkpoint
    else:
        dt = participant.get_max_time_step_size()
        U0 = participant.read_data(mesh_name, read_data_name, vertex_ids, dt)[0]
        t0 = t
        t = t0 + dt

# Stop coupling
participant.finalize()

import math

T_max = t
assert(math.isclose(T_max, 10))
I_analytical = lambda t: Io*np.cos(t*w0 + phi)
error = I0 - I_analytical(T_max)
print(f"error: {error}")