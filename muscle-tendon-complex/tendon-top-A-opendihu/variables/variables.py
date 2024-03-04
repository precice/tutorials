scenario_name = "tendon-top-A"

# time parameters
# ---------------
dt_elasticity = 1      # [ms] time step width for elasticity
end_time      = 20000   # [ms] simulation time
output_timestep_3D = 50  # [ms] output timestep

# setup
# -----
constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
force = 100.0           # [N] pulling force to the bottom 



# input files
# -----------

import os
opendihu_home = os.environ.get('OPENDIHU_HOME')
fiber_file = opendihu_home + "/examples/electrophysiology/input/left_biceps_brachii_tendon2a.bin"
cellml_file             = opendihu_home + "/examples/electrophysiology/input/2020_06_03_hodgkin-huxley_shorten_ocallaghan_davidson_soboleva_2007.cellml"
precice_config_file = "../precice-config.xml"
load_fiber_data = False             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.
debug_output = False                # verbose output in this python script, for debugging the domain decomposition
disable_firing_output = True        # Disables the initial list of fiber firings on the console to save some console space
paraview_output = False             # If the paraview output writer should be enabled
adios_output = False                # If the MegaMol/ADIOS output writer should be enabled
python_output = False               # If the Python output writer should be enabled
exfile_output = False               # If the Exfile output writer should be enabled# material parameters

# material parameters
# --------------------
tendon_material = "nonLinear"
rho = 10

# solvers
# -------
diffusion_solver_type = "cg"        # solver and preconditioner for the diffusion part of the Monodomain equation
diffusion_preconditioner_type = "none"      # preconditioner
potential_flow_solver_type = "gmres"        # solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_preconditioner_type = "none" # preconditioner
emg_solver_type = "cg"              # solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = False   #< If the initial guess for the emg linear system should be set to the previous solution

# partitioning
# ------------
# this has to match the total number of processes
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1

# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 50

mapping_tolerance = 0.1


# further internal variables that will be set by the helper.py script and used in the config in settings_fibers_emg.py
n_fibers_total = None
n_subdomains_xy = None
own_subdomain_coordinate_x = None
own_subdomain_coordinate_y = None
own_subdomain_coordinate_z = None
n_fibers_x = None
n_fibers_y = None
n_points_whole_fiber = None
n_points_3D_mesh_global_x = None
n_points_3D_mesh_global_y = None
n_points_3D_mesh_global_z = None
output_writer_fibers = None
output_writer_emg = None
output_writer_0D_states = None
states_output = False
parameters_used_as_algebraic = None
parameters_used_as_constant = None
parameters_initial_values = None
output_algebraic_index = None
output_state_index = None
nodal_stimulation_current = None
fiber_file_handle = None
fibers = None
fiber_distribution = None
firing_times = None
n_fibers_per_subdomain_x = None
n_fibers_per_subdomain_y = None
n_points_per_subdomain_z = None
z_point_index_start = None
z_point_index_end = None
meshes = None
potential_flow_dirichlet_bc = None
elasticity_dirichlet_bc = None
elasticity_neumann_bc = None
fibers_on_own_rank = None
n_fiber_nodes_on_subdomain = None
fiber_start_node_no = None
generate_linear_3d_mesh = False
generate_quadratic_3d_mesh = True
nx = None
ny = None
nz = None
constant_body_force = None
bottom_traction = None
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1
states_initial_values = []
enable_coupling = True
enable_force_length_relation = True
lambda_dot_scaling_factor = 1
mappings = None
vm_value_stimulated = None