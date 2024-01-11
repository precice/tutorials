from __future__ import division, print_function
import os
import sys
import argparse
import outputConfiguration as config
from thetaScheme import perform_partitioned_implicit_trapezoidal_rule_step, perform_partitioned_implicit_euler_step
import numpy as np
import tubePlotting
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from output import writeOutputToVTK
import precice

# physical properties of the tube
r0 = 1 / np.sqrt(np.pi)  # radius of the tube
a0 = r0**2 * np.pi  # cross sectional area
u0 = 10  # mean velocity
ampl = 3  # amplitude of varying velocity
frequency = 10  # frequency of variation
t_shift = 0  # temporal shift of variation
p0 = 0  # pressure at outlet
kappa = 100

L = 10  # length of tube/simulation domain
N = 100
dx = L / kappa
# helper function to create constant cross section


def velocity_in(t): return u0 + ampl * np.sin(frequency *
                                              (t + t_shift) * np.pi)  # inflow velocity


def crossSection0(N):
    return a0 * np.ones(N + 1)


parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml precice configuration file.",
                    nargs='?', type=str, default="../precice-config.xml")
parser.add_argument(
    "--enable-plot", help="Show a continuously updated plot of the tube while simulating.", action='store_true')
parser.add_argument("--write-video", help="Save a video of the simulation as 'writer_test.mp4'. \
                    NOTE: This requires 'enable_plot' to be active!", action='store_true')

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Did you forget adding the precice configuration file as an argument?")
    print("Try '$ python FluidSolver.py precice-config.xml'")
    quit()

plotting_mode = config.PlottingModes.VIDEO if args.enable_plot else config.PlottingModes.OFF
if args.write_video and not args.enable_plot:
    print("")
    print("To create a video it is required to enable plotting for this run.")
    print("Please supply both the '--enable-plot' and '--write-video' flags.")
    quit()
writeVideoToFile = True if args.write_video else False

print("Plotting Mode: {}".format(plotting_mode))

print("Starting Fluid Solver...")

print("N: " + str(N))

print("Configure preCICE...")
interface = precice.Participant("Fluid", args.configurationFileName, 0, 1)
print("preCICE configured...")

meshName = "Fluid-Nodes-Mesh"
crossSectionLengthName = "CrossSectionLength"
pressureName = "Pressure"

dimensions = interface.get_mesh_dimensions(meshName)

velocity = velocity_in(0) * np.ones(N + 1)
velocity_old = velocity_in(0) * np.ones(N + 1)
pressure = p0 * np.ones(N + 1)
pressure_old = p0 * np.ones(N + 1)
crossSectionLength = a0 * np.ones(N + 1)
crossSectionLength_old = a0 * np.ones(N + 1)

if plotting_mode == config.PlottingModes.VIDEO:
    fig, ax = plt.subplots(1)
    if writeVideoToFile:
        FFMpegWriter = manimation.writers['imagemagick']
        metadata = dict(title='PulseTube')
        writer = FFMpegWriter(fps=15, metadata=metadata)
        writer.setup(fig, "writer_test.mp4", 100)

grid = np.zeros([N + 1, dimensions])
grid[:, 0] = np.linspace(0, L, N + 1)  # x component
grid[:, 1] = 0  # y component, leave blank

vertexIDs = interface.set_mesh_vertices(meshName, grid)

if interface.requires_initial_data():
    interface.write_data(meshName, pressureName, vertexIDs, pressure)

t = 0


print("Fluid: init precice...")
# preCICE defines timestep size of solver via precice-config.xml
interface.initialize()

crossSectionLength = interface.read_data(
    meshName, crossSectionLengthName, vertexIDs, 0)

crossSectionLength_old = np.copy(crossSectionLength)
# initialize such that mass conservation is fulfilled
velocity_old = velocity_in(
    0) * crossSectionLength_old[0] * np.ones(N + 1) / crossSectionLength_old

print(crossSectionLength_old)

time_it = 0
while interface.is_coupling_ongoing():
    # When an implicit coupling scheme is used, checkpointing is required
    if interface.requires_writing_checkpoint():
        pass

    precice_dt = interface.get_max_time_step_size()

    crossSectionLength = interface.read_data(
        meshName, crossSectionLengthName, vertexIDs, precice_dt)
    velocity, pressure, success = perform_partitioned_implicit_euler_step(
        velocity_old, pressure_old, crossSectionLength_old, crossSectionLength, dx, precice_dt, velocity_in(
            t + precice_dt), custom_coupling=True)
    interface.write_data(meshName, pressureName, vertexIDs, pressure)
    interface.advance(precice_dt)

    # i.e. not yet converged
    if interface.requires_reading_checkpoint():
        pass
    else:  # converged, timestep complete
        t += precice_dt
        if plotting_mode is config.PlottingModes.VIDEO:
            tubePlotting.doPlotting(
                ax, crossSectionLength_old, velocity_old, pressure_old, dx, t)
            if writeVideoToFile:
                writer.grab_frame()
            ax.cla()
        velocity_old = np.copy(velocity)
        pressure_old = np.copy(pressure)
        crossSectionLength_old = np.copy(crossSectionLength)
        writeOutputToVTK(time_it, "out_fluid_", dx, datanames=["velocity", "pressure", "diameter"], data=[
            velocity_old, pressure_old, crossSectionLength_old])
        time_it += 1

print("Exiting FluidSolver")

if plotting_mode is config.PlottingModes.VIDEO and writeVideoToFile:
    writer.finish()
