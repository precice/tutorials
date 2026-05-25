from mesa.examples.advanced.wolf_sheep.agents import GrassPatch
from mesa.examples.advanced.wolf_sheep.model import WolfSheep
from mesa.examples.advanced.wolf_sheep.model import WolfSheepScenario
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import precice

min_depth_for_grass = 0.2

# Use dx of Raster Model Grid from Soil Creep for solver_dt calculation
dx = 1.0
fast_creep = 0.1

# Initialize Wolf-Sheep-Grass model
wss = WolfSheepScenario(
    initial_sheep=20,
    initial_wolves=10,
    grass=True,
    grass_regrowth_time=15  # give grass a fighting chance...
)
ws = WolfSheep(wss)

ground_cover_cmap = copy.copy(mpl.colormaps["YlGn"])


def generate_grass_map(model):
    grass_map = np.zeros((model.grid.width, model.grid.height))
    for cell in model.grid:
        (x, y) = cell.coordinate
        cell_content = cell.agents
        for agent in cell_content:
            if type(agent) is GrassPatch:
                if agent.fully_grown:
                    grass_map[x][y] = 2
                else:
                    grass_map[x][y] = 1
    return grass_map


def plot_grass_map(grass_map):
    plt.imshow(grass_map, interpolation="nearest", cmap=ground_cover_cmap)
    plt.colorbar()


def limit_grass_by_soil(wsg_model, soil, min_soil_depth):
    soilmatrix = soil.reshape((wsg_model.width, wsg_model.height))
    for cell in wsg_model.grid:
        (x, y) = cell.coordinate
        cell_content = cell.agents
        if soilmatrix[x][y] < min_soil_depth:
            for agent in cell_content:
                if type(agent) is GrassPatch:
                    agent.fully_grown = False


width = ws.grid.width
height = ws.grid.height

# preCICE setup
participant_name = "Wolf-Sheep-Grass"
config_file_name = "../precice-config.xml"
solver_process_index = 0
solver_process_size = 1
participant = precice.Participant(participant_name, config_file_name, solver_process_index, solver_process_size)

positions = [[x, y] for x in range(width) for y in range(height)]
vertex_gm_ids = participant.set_mesh_vertices("Wolf-Sheep-Grass-Mesh", positions)
vertex_soil_ids = participant.set_mesh_vertices("Soil-Grass-Mesh", positions)

soil = np.zeros([width * height])

if participant.requires_initial_data():
    gm = generate_grass_map(ws)
    participant.write_data("Wolf-Sheep-Grass-Mesh", "Grass", vertex_gm_ids, gm.flatten())

participant.initialize()

while participant.is_coupling_ongoing():
    solver_dt = 0.2 * dx * dx / fast_creep
    precice_dt = participant.get_max_time_step_size()
    dt = np.minimum(solver_dt, precice_dt)

    soil = participant.read_data("Soil-Grass-Mesh", "Soil", vertex_soil_ids, dt)

    # Update the grass cover
    limit_grass_by_soil(ws, soil, min_depth_for_grass)

    # Run the W-S-G model
    ws.step()

    gm = generate_grass_map(ws)
    participant.write_data("Wolf-Sheep-Grass-Mesh", "Grass", vertex_gm_ids, gm.flatten())

    participant.advance(dt)

participant.finalize()
