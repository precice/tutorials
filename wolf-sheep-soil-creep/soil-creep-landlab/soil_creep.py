import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from landlab import RasterModelGrid, imshow_grid
from landlab.components import LinearDiffuser
import precice

initial_soil_depth = 0.3
hstar = 0.2
fast_creep = 0.1
slow_creep = 0.001

ground_cover_cmap = copy.copy(mpl.colormaps["YlGn"])

# Create a grid the same size as the W-S-G model's grid
width = 20
height = 20

rmg = RasterModelGrid((width, height))

# Create elevation field and have it slope down to the south at 10% gradient
elev = rmg.add_zeros("topographic__elevation", at="node")
elev[:] = 0.1 * rmg.y_of_node

# Have one open boundary on the south side
rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)

# Remember the starting elevation so we can calculate cumulative erosion/deposition
initial_elev = np.zeros(rmg.number_of_nodes)
initial_elev[:] = elev

# Also remember the elevation of the prior time step, so we can difference
prior_elev = np.zeros(rmg.number_of_nodes)

# Create a field for the creep coefficient, and set parameters for two
# rates: slow (full grass cover) and fast (partial or "eaten" grass cover)
creep_coef = rmg.add_zeros("creep_coefficient", at="node")

# Create a soil-thickness field
soil = rmg.add_zeros("soil__depth", at="node")
soil[:] = initial_soil_depth

# Instantiate a LinearDiffuser (soil creep) Landlab component
diffuser = LinearDiffuser(rmg, linear_diffusivity=creep_coef)

# preCICE setup
participant_name = "Soil-Creep"
config_file_name = "../precice-config.xml"
solver_process_index = 0
solver_process_size = 1
participant = precice.Participant(participant_name, config_file_name, solver_process_index, solver_process_size)

positions = [[x, y] for x in range(width) for y in range(height)]
vertex_gm_ids = participant.set_mesh_vertices("Soil-Creep-Mesh", positions)
vertex_soil_ids = participant.set_mesh_vertices("Soil-Depth-Mesh", positions)

participant.initialize()

while participant.is_coupling_ongoing():
    solver_dt = 0.2 * rmg.dx * rmg.dx / fast_creep
    precice_dt = participant.get_max_time_step_size()
    dt = np.minimum(solver_dt, precice_dt)

    gm = participant.read_data("Soil-Creep-Mesh", "Grass", vertex_gm_ids, dt)

    # Assign the higher creep coefficient to cells where the grass has
    # been eaten and not yet recovered; the slower value is assigned to
    # "fully grown" grass patches.
    creep_coef[gm.flatten() == 1] = fast_creep
    creep_coef[gm.flatten() == 2] = slow_creep

    # Adjust the creep coefficient to account for soil depth
    creep_coef *= 1.0 - np.exp(-soil / hstar)

    # Run the soil-creep model
    prior_elev[:] = elev
    diffuser.run_one_step(dt)

    # Update the soil cover
    soil += elev - prior_elev

    participant.write_data("Soil-Depth-Mesh", "Soil", vertex_soil_ids, soil)

    participant.advance(dt)

participant.finalize()

# Calculate and plot the erosion/deposition patterns
ero_dep = elev - initial_elev
maxchange = np.amax(np.abs(ero_dep))
imshow_grid(
    rmg,
    ero_dep,
    vmin=-maxchange,
    vmax=maxchange,
    cmap="coolwarm_r",
    colorbar_label="Depth of soil accumulation (+) or loss (-), m",
)
plt.show()

# Soil thickness
imshow_grid(rmg, soil, colorbar_label="Soil thickness, m")
plt.show()

# Ground cover
imshow_grid(
    rmg, gm, cmap=ground_cover_cmap, colorbar_label="Ground cover (1 = bare, 2 = grass)"
)
plt.show()
