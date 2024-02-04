from __future__ import division

import argparse
import numpy as np
from numpy.linalg import eig
import precice
from enum import Enum
import csv
import os

import problemDefinition


class Scheme(Enum):
    NEWMARK_BETA = "Newmark_beta"
    GENERALIZED_ALPHA = "generalized_alpha"


class Participant(Enum):
    MASS_LEFT = "Mass-Left"
    MASS_RIGHT = "Mass-Right"


parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str, choices=[p.value for p in Participant])
parser.add_argument("-ts", "--time-stepping", help="Time stepping scheme being used.", type=str,
                    choices=[s.value for s in Scheme], default=Scheme.NEWMARK_BETA.value)
args = parser.parse_args()

participant_name = args.participantName

if participant_name == Participant.MASS_LEFT.value:
    write_data_name = 'Displacement-Left'
    read_data_name = 'Displacement-Right'
    mesh_name = 'Mass-Left-Mesh'

    this_mass = problemDefinition.MassLeft
    this_spring = problemDefinition.SpringLeft
    connecting_spring = problemDefinition.SpringMiddle
    other_mass = problemDefinition.MassRight

elif participant_name == Participant.MASS_RIGHT.value:
    read_data_name = 'Displacement-Left'
    write_data_name = 'Displacement-Right'
    mesh_name = 'Mass-Right-Mesh'

    this_mass = problemDefinition.MassRight
    this_spring = problemDefinition.SpringRight
    connecting_spring = problemDefinition.SpringMiddle
    other_mass = problemDefinition.MassLeft

else:
    raise Exception(f"wrong participant name: {participant_name}")

mass = this_mass.m
stiffness = this_spring.k + connecting_spring.k
u0, v0, f0, d_dt_f0 = this_mass.u0, this_mass.v0, connecting_spring.k * other_mass.u0, connecting_spring.k * other_mass.v0

num_vertices = 1  # Number of vertices

solver_process_index = 0
solver_process_size = 1

configuration_file_name = "../precice-config.xml"

participant = precice.Participant(participant_name, configuration_file_name, solver_process_index, solver_process_size)

dimensions = participant.get_mesh_dimensions(mesh_name)

vertex = np.zeros(dimensions)
read_data = np.zeros(num_vertices)
write_data = connecting_spring.k * u0 * np.ones(num_vertices)

vertex_ids = [participant.set_mesh_vertex(mesh_name, vertex)]

if participant.requires_initial_data():
    participant.write_data(mesh_name, write_data_name, vertex_ids, write_data)

participant.initialize()
precice_dt = participant.get_max_time_step_size()
my_dt = precice_dt  # use my_dt < precice_dt for subcycling

# Initial Conditions
a0 = (f0 - stiffness * u0) / mass
u = u0
v = v0
a = a0
t = 0

# Generalized Alpha Parameters
if args.time_stepping == Scheme.GENERALIZED_ALPHA.value:
    alpha_f = 0.4
    alpha_m = 0.2
elif args.time_stepping == Scheme.NEWMARK_BETA.value:
    alpha_f = 0.0
    alpha_m = 0.0
m = 3 * [None]  # will be computed for each timestep depending on dt
gamma = 0.5 - alpha_m + alpha_f
beta = 0.25 * (gamma + 0.5)

positions = []
velocities = []
times = []

u_write = [u]
v_write = [v]
t_write = [t]

while participant.is_coupling_ongoing():
    if participant.requires_writing_checkpoint():
        u_cp = u
        v_cp = v
        a_cp = a
        t_cp = t
        # store data for plotting and postprocessing
        positions += u_write
        velocities += v_write
        times += t_write

    # compute time step size for this time step
    precice_dt = participant.get_max_time_step_size()
    dt = np.min([precice_dt, my_dt])
    read_time = (1 - alpha_f) * dt
    read_data = participant.read_data(mesh_name, read_data_name, vertex_ids, read_time)
    f = connecting_spring.k * read_data[0]

    # do generalized alpha step
    m[0] = (1 - alpha_m) / (beta * dt**2)
    m[1] = (1 - alpha_m) / (beta * dt)
    m[2] = (1 - alpha_m - 2 * beta) / (2 * beta)
    k_bar = stiffness * (1 - alpha_f) + m[0] * mass
    u_new = (f - alpha_f * stiffness * u + mass * (m[0] * u + m[1] * v + m[2] * a)) / k_bar
    a_new = 1.0 / (beta * dt**2) * (u_new - u - dt * v) - (1 - 2 * beta) / (2 * beta) * a
    v_new = v + dt * ((1 - gamma) * a + gamma * a_new)
    t_new = t + dt

    write_data = [u_new]

    participant.write_data(mesh_name, write_data_name, vertex_ids, write_data)

    participant.advance(dt)

    if participant.requires_reading_checkpoint():
        u = u_cp
        v = v_cp
        a = a_cp
        t = t_cp
        # empty buffers for next window
        u_write = []
        v_write = []
        t_write = []

    else:
        u = u_new
        v = v_new
        a = a_new
        t = t_new

        # write data to buffers
        u_write.append(u)
        v_write.append(v)
        t_write.append(t)

# store final result
u = u_new
v = v_new
a = a_new
u_write.append(u)
v_write.append(v)
t_write.append(t)
positions += u_write
velocities += v_write
times += t_write

participant.finalize()

# print errors
error = np.max(abs(this_mass.u_analytical(np.array(times)) - np.array(positions)))
print("Error w.r.t analytical solution:")
print(f"{my_dt},{error}")

# output trajectory
if not os.path.exists("output"):
    os.makedirs("output")

with open(f'output/trajectory-{participant_name}.csv', 'w') as file:
    csv_write = csv.writer(file, delimiter=';')
    csv_write.writerow(['time', 'position', 'velocity'])
    for t, u, v in zip(times, positions, velocities):
        csv_write.writerow([t, u, v])
