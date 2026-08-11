import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import precice

class LandlabLinearDiffuserClone:
    CORE = 0 # regular node
    FIXED_VALUE = 1 # Dirichlet fixed elevation
    CLOSED = 4 # Neumann zero-flux

    ALPHA = 0.15  # time-step stability factor

    def __init__(self, z, D, dx):
        self.z = z # grid
        self.D = D # linear diffusivity
        self.dx = dx # node spacing
        self.ny, self.nx = z.shape # grid height and width
        self.id = np.arange(self.nx * self.ny).reshape(self.ny, self.nx)

        self.status = np.full((self.ny, self.nx), self.FIXED_VALUE, dtype=np.uint8)
        self.status[1:-1, 1:-1] = self.CORE

        # Have one open boundary on the south side
        # Corresponds to Landlab's rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        # with the boundary order: right, top, left, bottom
        self.status[:, -1] = self.CLOSED
        self.status[-1, :] = self.CLOSED
        self.status[:, 0] = self.CLOSED

        self._build_links()

    def _build_links(self): # Initalize list of active links: between two core nodes or between core and fixed
        ny, nx = self.ny, self.nx
        links = []

        for j in range(ny):
            for i in range(nx - 1):
                links.append((self.id[j, i], self.id[j, i + 1]))

        for j in range(ny - 1):
            for i in range(nx):
                links.append((self.id[j, i], self.id[j + 1, i]))

        status = self.status.ravel()
        active_links = []

        for n1, n2 in links:
            s1 = status[n1]
            s2 = status[n2]

            if s1 == self.CLOSED or s2 == self.CLOSED:
                continue

            if s1 != self.CORE and s2 != self.CORE:
                continue

            active_links.append((n1, n2))

        self.links = np.array(active_links, dtype=int)

    def step(self, dt):
        D = self.D.ravel()

        active_D = np.array([
            np.maximum(D[n1], D[n2])
            for n1, n2 in self.links
        ])

        internal_dt = ALPHA * self.dx * self.dx / np.nanmax(active_D)

        repeats = int(dt // internal_dt)
        remainder = dt - repeats * internal_dt

        for _ in range(repeats):
            self._substep(internal_dt)

        if remainder > 0:
            self._substep(remainder)

    def _substep(self, dt):
        z = self.z.ravel()
        D = self.D.ravel()
        status = self.status.ravel()

        dzdt = np.zeros_like(z)

        for n1, n2 in self.links: # Iterate over all active links (tuples)
            Dk = np.maximum(D[n1], D[n2]) # max link diffusivity
            q = -Dk * (z[n2] - z[n1]) / self.dx

            if status[n1] == self.CORE:
                dzdt[n1] -= q / self.dx
            if status[n2] == self.CORE:
                dzdt[n2] += q / self.dx

        z[status == self.CORE] += dt * dzdt[status == self.CORE]
        self.z = z.reshape(self.ny, self.nx)

initial_soil_depth = 0.3

# Set the soil-thickness scale for limiting creep where little soil is available
hstar = 0.2

# Set parameters for two soil creep coefficients: slow (full grass cover) and fast (partial or "eaten" grass cover)
fast_creep = 0.1
slow_creep = 0.001

# Set grid size (same as W-S-G model's grid)
# Adjust refinement by changing number of nodes ny/nx (e.g. double refinement by setting ny = nx = 40).
ny = nx = 20
length_x = 19
length_y = 19
dx = length_x/(nx-1)

solver_dt = 0.2 * dx * dx / fast_creep

# Grass map field for plotting
gm = np.zeros((ny,nx))

# Create elevation field and have it slope down to the south at 10% gradient
y = np.arange(ny)[:, None] * np.ones((ny, nx)) * dx
elev = 0.1 * y

# Remember the starting elevation so we can calculate cumulative erosion/deposition
initial_elev = np.zeros((ny,nx))
initial_elev[:] = elev

# Also remember the elevation of the prior time step, so we can difference
prior_elev = np.zeros((ny,nx))

# Create a field for the creep coefficient
creep_coef = np.zeros((ny,nx))

# Create a soil-thickness field
soil = np.full((ny, nx), initial_soil_depth)

# Instantiate a LinearDiffuser (soil creep) Landlab component
diffuser = LandlabLinearDiffuserClone(
    z=elev,
    D=creep_coef,
    dx=dx,
)

# preCICE setup
participant_name = "Soil-Creep"
config_file_name = "../precice-config.xml"
solver_process_index = 0
solver_process_size = 1
participant = precice.Participant(participant_name, config_file_name, solver_process_index, solver_process_size)

x = np.linspace(0, length_x, nx)
y = np.linspace(0, length_y, ny)

positions = [[x0, y0] for y0 in y for x0 in x]
vertex_ids = participant.set_mesh_vertices("Soil-Creep-Mesh", positions)

participant.initialize()

while participant.is_coupling_ongoing():
    precice_dt = participant.get_max_time_step_size()
    dt = np.minimum(solver_dt, precice_dt)

    gm_flat = participant.read_data("Soil-Creep-Mesh", "Grass", vertex_ids, dt)
    gm[:,:] = np.asarray(gm_flat).reshape(ny, nx)

    # Assign the higher creep coefficient to cells where the grass has
    # been eaten and not yet recovered; the slower value is assigned to
    # "fully grown" grass patches.
    creep_coef[gm == 1] = fast_creep
    creep_coef[gm == 2] = slow_creep

    # Adjust the creep coefficient to account for soil depth
    creep_coef[:, :] *= 1.0 - np.exp(-soil / hstar)

    # Remember the current elevation
    prior_elev[:, :] = elev

    # Run the soil-creep model
    diffuser.step(dt)

    # Update elevation
    elev = diffuser.z

    # Update the soil cover
    soil += elev - prior_elev

    participant.write_data("Soil-Creep-Mesh", "Soil", vertex_ids, soil.ravel())

    participant.advance(dt)

participant.finalize()


# Mask for closed boundaries (mimic Landlab)
closed = np.zeros((ny, nx), dtype=bool)

closed[1:-1, 0] = True
closed[1:-1, -1] = True
closed[-1, :] = True

# Calculate and plot the erosion/deposition patterns
plt.figure()
ero_dep = elev - initial_elev
maxchange = np.max(np.abs(ero_dep))

cmap = mpl.colormaps["coolwarm_r"].copy()
cmap.set_bad("black")

plt.figure()
plt.imshow(
    np.ma.array(ero_dep, mask=closed),
    origin="lower",
    vmin=-maxchange,
    vmax=maxchange,
    cmap=cmap
)
plt.colorbar(label="Depth of soil accumulation (+) or loss (-), m")
plt.savefig("output/erosion_deposition_patterns.png")
plt.close()

# Soil thickness
cmap = mpl.colormaps["pink"].copy()
cmap.set_bad("black")

plt.figure()
plt.imshow(
    np.ma.array(soil, mask=closed),
    origin="lower",
    cmap=cmap
)
plt.colorbar(label="Soil thickness, m")
plt.savefig("output/soil_thickness.png")
plt.close()

# Ground cover
cmap = mpl.colormaps["YlGn"].copy()
cmap.set_bad("black")

plt.figure()
plt.imshow(
    np.ma.array(gm, mask=closed),
    origin="lower",
    cmap=cmap
)
plt.colorbar(label="Ground cover (1 = bare, 2 = grass)")
plt.savefig("output/grass_map.png")
plt.close()
