# material parameters
# --------------------

c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b = 1.075e-2               # [N/cm^2] anisotropy parameter
d = 9.1733                 # [-] anisotropy parameter
material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                          # maximum stress [N/cm^2]
Conductivity = 3.828                # sigma, conductivity [mS/cm]
Am = 500.0                          # surface area to volume ratio [cm^-1]
Cm = 0.58                           # membrane capacitance [uF/cm^2]
# not used [cm], this will later be used to specify a variance of positions of the innervation point at the fibers
innervation_zone_width = 0.
rho = 10

# solvers
# -------
diffusion_solver_type = "cg"        # solver and preconditioner for the diffusion part of the Monodomain equation
diffusion_preconditioner_type = "none"      # preconditioner
# solver and preconditioner for an initial Laplace flow on the domain, from which fiber directions are determined
potential_flow_solver_type = "gmres"
potential_flow_preconditioner_type = "none"  # preconditioner
# solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_solver_type = "cg"
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = False  # < If the initial guess for the emg linear system should be set to the previous solution

# timing parameters
# -----------------
end_time = 1000.0                   # [ms] end time of the simulation
# [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency = 100 * 1e-3
dt_0D = 1e-3                        # [ms] timestep width of ODEs
dt_1D = 1.5e-3                      # [ms] timestep width of diffusion
dt_splitting = 3e-3                 # [ms] overall timestep width of strang splitting
# [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
dt_3D = 1e0
output_timestep = 1e0               # [ms] timestep for output files
output_timestep_3D_emg = 1e0        # [ms] timestep for output files
output_timestep_3D = 1e0            # [ms] timestep for output files
output_timestep_fibers = 1e0        # [ms] timestep for output files
activation_start_time = 0           # [ms] time when to start checking for stimulation

# input files
# -----------

import os
input_dir = os.environ.get('OPENDIHU_INPUT_DIR')
fiber_file = input_dir + "/left_biceps_brachii_9x9fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
# use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
firing_times_file = input_dir + "/MU_firing_times_always.txt"
fiber_distribution_file = input_dir + "/MU_fibre_distribution_10MUs.txt"
cellml_file = input_dir + "/2020_06_03_hodgkin-huxley_shorten_ocallaghan_davidson_soboleva_2007.cellml"
firing_times_file = input_dir + "/MU_firing_times_real.txt"
precice_config_file = "../precice-config.xml"
# If the fiber geometry data should be loaded completely in the python
# script. If True, this reads the binary file and assigns the node
# positions in the config. If False, the C++ code will read the binary
# file and only extract the local node positions. This is more performant
# for highly parallel runs.
load_fiber_data = False
debug_output = False                # verbose output in this python script, for debugging the domain decomposition
disable_firing_output = True        # Disables the initial list of fiber firings on the console to save some console space
paraview_output = False             # If the paraview output writer should be enabled
adios_output = False                # If the MegaMol/ADIOS output writer should be enabled
python_output = False               # If the Python output writer should be enabled
exfile_output = False               # If the Exfile output writer should be enabled

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

# scenario name for log file
scenario_name = ""

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
# These functions can be redefined differently in a custom variables script


def get_am(fiber_no, mu_no):
    return Am


def get_cm(fiber_no, mu_no):
    return Cm


def get_conductivity(fiber_no, mu_no):
    return Conductivity


def get_specific_states_call_frequency(fiber_no, mu_no):
    return stimulation_frequency


def get_specific_states_frequency_jitter(fiber_no, mu_no):
    return [0]


def get_specific_states_call_enable_begin(fiber_no, mu_no):
    return activation_start_time


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
