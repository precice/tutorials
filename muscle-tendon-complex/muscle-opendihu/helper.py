# Multiple 1D fibers (monodomain) with 3D contraction, biceps geometry
# This is a helper script that sets a lot of the internal variables which are all defined in variables.py
#
# if variables.fiber_file=cuboid.bin, it uses a small cuboid test example

import numpy as np
import pickle
import sys
import os
import struct
import argparse
# sys.path.insert(0, '..')
import variables    # file variables.py
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# generate cuboid fiber file
if "cuboid.bin" in variables.fiber_file:

    if variables.n_fibers_y is None:
        variables.n_fibers_x = 4
        variables.n_fibers_y = variables.n_fibers_x
        variables.n_points_whole_fiber = 20

    size_x = variables.n_fibers_x * 0.1
    size_y = variables.n_fibers_y * 0.1
    size_z = variables.n_points_whole_fiber / 100.

    if rank_no == 0:
        print(
            "create cuboid.bin with size [{},{},{}], n points [{},{},{}]".format(
                size_x,
                size_y,
                size_z,
                variables.n_fibers_x,
                variables.n_fibers_y,
                variables.n_points_whole_fiber))

        # write header
        with open(variables.fiber_file, "wb") as outfile:

            # write header
            header_str = "opendihu self-generated cuboid  "
            outfile.write(struct.pack('32s', bytes(header_str, 'utf-8')))   # 32 bytes
            outfile.write(struct.pack('i', 40))  # header length
            outfile.write(struct.pack('i', variables.n_fibers_x * variables.n_fibers_y))   # n_fibers
            outfile.write(struct.pack('i', variables.n_points_whole_fiber))   # variables.n_points_whole_fiber
            outfile.write(struct.pack('i', 0))   # nBoundaryPointsXNew
            outfile.write(struct.pack('i', 0))   # nBoundaryPointsZNew
            outfile.write(struct.pack('i', 0))   # nFineGridFibers_
            outfile.write(struct.pack('i', 1))   # nRanks
            outfile.write(struct.pack('i', 1))   # nRanksZ
            outfile.write(struct.pack('i', 0))   # nFibersPerRank
            outfile.write(struct.pack('i', 0))   # date

            # loop over points
            for y in range(variables.n_fibers_y):
                for x in range(variables.n_fibers_x):
                    for z in range(variables.n_points_whole_fiber):
                        point = [x * (float)(size_x) / (variables.n_fibers_x), y * (float)(size_y) /
                                 (variables.n_fibers_y), z * (float)(size_z) / (variables.n_points_whole_fiber)]
                        outfile.write(struct.pack('3d', point[0], point[1], point[2]))   # data point

# output diffusion solver type
if rank_no == 0:
    print("diffusion solver type: {}".format(variables.diffusion_solver_type))

variables.load_fiber_data = False   # load all local node positions from fiber_file, in order to infer partitioning for fat_layer mesh

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x,
    variables.n_subdomains_y,
    variables.n_subdomains_z,
    variables.fiber_file,
    variables.load_fiber_data,
    variables.sampling_stride_x,
    variables.sampling_stride_y,
    variables.sampling_stride_z,
    variables.generate_linear_3d_mesh,
    variables.generate_quadratic_3d_mesh)
[variables.meshes,
    variables.own_subdomain_coordinate_x,
    variables.own_subdomain_coordinate_y,
    variables.own_subdomain_coordinate_z,
    variables.n_fibers_x,
    variables.n_fibers_y,
 variables.n_points_whole_fiber] = result

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

# create mappings between meshes
# variables.mappings_between_meshes = {"MeshFiber_{}".format(i) : "3Dmesh" for i in range(variables.n_fibers_total)}
variables.mappings_between_meshes = {"MeshFiber_{}".format(
    i): {"name": "3Dmesh", "xiTolerance": 1e-3} for i in range(variables.n_fibers_total)}

# a higher tolerance includes more fiber dofs that may be almost out of the 3D mesh
variables.mappings_between_meshes = {
    "MeshFiber_{}".format(i): {
        "name": "3Dmesh_quadratic",
        "xiTolerance": variables.mapping_tolerance,
        "enableWarnings": False,
        "compositeUseOnlyInitializedMappings": False,
        "fixUnmappedDofs": True,
        "defaultValue": 0,
    } for i in range(variables.n_fibers_total)
}
# set output writer
variables.output_writer_fibers = []
variables.output_writer_elasticity = []
variables.output_writer_emg = []
variables.output_writer_0D_states = []

subfolder = ""
if variables.paraview_output:
    if variables.adios_output:
        subfolder = "paraview/"
    variables.output_writer_emg.append({"format": "Paraview",
                                        "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D_emg),
                                        "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg",
                                        "binary": True,
                                        "fixedFormat": False,
                                        "combineFiles": True,
                                        "fileNumbering": "incremental"})
    variables.output_writer_elasticity.append({"format": "Paraview",
                                               "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D),
                                               "filename": "out/" + subfolder + variables.scenario_name + "/elasticity",
                                               "binary": True,
                                               "fixedFormat": False,
                                               "combineFiles": True,
                                               "fileNumbering": "incremental"})
    variables.output_writer_fibers.append({"format": "Paraview",
                                           "outputInterval": int(1. / variables.dt_splitting * variables.output_timestep_fibers),
                                           "filename": "out/" + subfolder + variables.scenario_name + "/fibers",
                                           "binary": True,
                                           "fixedFormat": False,
                                           "combineFiles": True,
                                           "fileNumbering": "incremental"})
    if variables.states_output:
        variables.output_writer_0D_states.append({"format": "Paraview",
                                                  "outputInterval": 1,
                                                  "filename": "out/" + subfolder + variables.scenario_name + "/0D_states",
                                                  "binary": True,
                                                  "fixedFormat": False,
                                                  "combineFiles": True,
                                                  "fileNumbering": "incremental"})

if variables.adios_output:
    if variables.paraview_output:
        subfolder = "adios/"
    variables.output_writer_emg.append({"format": "MegaMol",
                                        "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D_emg),
                                        "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg",
                                        "useFrontBackBuffer": False,
                                        "combineNInstances": 1,
                                        "fileNumbering": "incremental"})
    variables.output_writer_elasticity.append({"format": "MegaMol",
                                               "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D),
                                               "filename": "out/" + subfolder + variables.scenario_name + "/elasticity",
                                               "useFrontBackBuffer": False,
                                               "fileNumbering": "incremental"})
    variables.output_writer_fibers.append({"format": "MegaMol",
                                           "outputInterval": int(1. / variables.dt_splitting * variables.output_timestep_fibers),
                                           "filename": "out/" + subfolder + variables.scenario_name + "/fibers",
                                           "combineNInstances": variables.n_subdomains_xy,
                                           "useFrontBackBuffer": False,
                                           "fileNumbering": "incremental"})
    # variables.output_writer_fibers.append({"format": "MegaMol",
    # "outputInterval":
    # int(1./variables.dt_splitting*variables.output_timestep_fibers),
    # "filename": "out/" + variables.scenario_name + "/fibers",
    # "combineNInstances": 1, "useFrontBackBuffer": False, "fileNumbering":
    # "incremental"}

if variables.python_output:
    if variables.adios_output:
        subfolder = "python/"
    variables.output_writer_emg.append({"format": "PythonFile",
                                        "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D_emg),
                                        "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg",
                                        "binary": True,
                                        "fileNumbering": "incremental"})
    variables.output_writer_elasticity.append({"format": "PythonFile",
                                               "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D),
                                               "filename": "out/" + subfolder + variables.scenario_name + "/elasticity",
                                               "binary": True,
                                               "fileNumbering": "incremental"})
    variables.output_writer_fibers.append({"format": "PythonFile",
                                           "outputInterval": int(1. / variables.dt_splitting * variables.output_timestep_fibers),
                                           "filename": "out/" + subfolder + variables.scenario_name + "/fibers",
                                           "binary": True,
                                           "fileNumbering": "incremental"})

if variables.exfile_output:
    if variables.adios_output:
        subfolder = "exfile/"
    variables.output_writer_emg.append({"format": "Exfile",
                                        "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D_emg),
                                        "filename": "out/" + subfolder + variables.scenario_name + "/hd_emg",
                                        "fileNumbering": "incremental"})
    variables.output_writer_elasticity.append({"format": "Exfile",
                                               "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D),
                                               "filename": "out/" + subfolder + variables.scenario_name + "/elasticity",
                                               "fileNumbering": "incremental"})
    variables.output_writer_fibers.append({"format": "Exfile",
                                           "outputInterval": int(1. / variables.dt_splitting * variables.output_timestep_fibers),
                                           "filename": "out/" + subfolder + variables.scenario_name + "/fibers",
                                           "fileNumbering": "incremental"})

# set variable mappings for cellml model
if "hodgkin-huxley_shorten_ocallaghan_davidson_soboleva_2007.cellml" in variables.cellml_file:   # hodgkin huxley membrane with fatigue from shorten
    # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
    variables.mappings = {
        ("parameter", 0): ("constant", "membrane/i_Stim"),      # parameter 0 is I_stim
        ("parameter", 1): ("constant", "razumova/l_hs"),        # parameter 1 is fiber stretch λ
        ("parameter", 2): ("constant", "razumova/velocity"),    # parameter 2 is fiber contraction velocity \dot{λ}
        ("connectorSlot", "vm"): ("state", "membrane/V"),           # expose state Vm to the operator splitting
        ("connectorSlot", "stress"): ("algebraic", "razumova/activestress"),  # expose algebraic γ to the operator splitting
        ("connectorSlot", "alpha"): ("algebraic", "razumova/activation"),  # expose algebraic α to the operator splitting

        # connect output "lamda" of mechanics solver to parameter 1 (l_hs)
        ("connectorSlot", "lambda"): "razumova/l_hs",
        # connect output "ldot"  of mechanics solver to parameter 2 (rel_velo)
        ("connectorSlot", "ldot"): "razumova/velocity",
    }
    # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
    variables.parameters_initial_values = [0, 1, 0]
    variables.nodal_stimulation_current = 40.                           # not used
    # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
    variables.vm_value_stimulated = 40.
    # disable computation of force-length relation in opendihu, as it is carried out in CellML model
    variables.enable_force_length_relation = False
    # scaling factor to convert dimensionless contraction velocity to
    # shortening velocity, velocity = factor*\dot{lambda}
    variables.lambda_dot_scaling_factor = 7.815e-05

elif "hodgkin_huxley" in variables.cellml_file:
    # parameters: I_stim
    variables.mappings = {
        ("parameter", 0): ("constant", "membrane/i_Stim"),     # parameter 0 is constant 2 = I_stim
        ("connectorSlot", "vm"): "membrane/V",                        # expose state 0 = Vm to the operator splitting
    }
    variables.parameters_initial_values = [0.0]                         # initial value for stimulation current
    variables.nodal_stimulation_current = 40.                           # not used
    # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
    variables.vm_value_stimulated = 20.

elif "shorten" in variables.cellml_file:
    # parameters: stimulation current I_stim, fiber stretch λ
    variables.mappings = {
        ("parameter", 0): ("algebraic", "wal_environment/I_HH"),  # parameter is algebraic 32
        # parameter is constant 65, fiber stretch λ, this indicates how much the
        # fiber has stretched, 1 means no extension
        ("parameter", 1): ("constant", "razumova/L_x"),
        ("connectorSlot", "vm"): "wal_environment/vS",                 # expose state 0 = Vm to the operator splitting
    }
    variables.parameters_initial_values = [0.0, 1.0]                    # stimulation current I_stim, fiber stretch λ
    variables.nodal_stimulation_current = 1200.                         # not used
    # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
    variables.vm_value_stimulated = 40.

elif "slow_TK_2014" in variables.cellml_file:   # this is (3a, "MultiPhysStrain", old tomo mechanics) in OpenCMISS
    # parameters: I_stim, fiber stretch λ
    variables.mappings = {
        ("parameter", 0): ("constant", "wal_environment/I_HH"),  # parameter 0 is constant 54 = I_stim
        ("parameter", 1): ("constant", "razumova/L_S"),         # parameter 1 is constant 67 = fiber stretch λ
        ("connectorSlot", "vm"): "wal_environment/vS",                 # expose state 0 = Vm to the operator splitting
        # expose algebraic 12 = γ to the operator splitting
        ("connectorSlot", "stress"): "razumova/stress",
    }
    # wal_environment/I_HH = I_stim, razumova/L_S = λ
    variables.parameters_initial_values = [0.0, 1.0]
    variables.nodal_stimulation_current = 40.                           # not used
    # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
    variables.vm_value_stimulated = 40.

# this is (3, "MultiPhysStrain", numerically more stable) in OpenCMISS, this only computes A1,A2,x1,x2 not the stress
elif "Aliev_Panfilov_Razumova_2016_08_22" in variables.cellml_file:
    # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
    variables.mappings = {
        ("parameter", 0): ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
        ("parameter", 1): ("constant", "Razumova/l_hs"),        # parameter 1 is constant 8 = fiber stretch λ
        # parameter 2 is constant 9 = fiber contraction velocity \dot{λ}
        ("parameter", 2): ("constant", "Razumova/velo"),
        ("connectorSlot", "vm"): "Aliev_Panfilov/V_m",                 # expose state 0 = Vm to the operator splitting
        # expose algebraic 0 = γ to the operator splitting
        ("connectorSlot", "stress"): "Razumova/sigma",
    }
    # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/velo = \dot{λ}
    variables.parameters_initial_values = [0, 1, 0]
    variables.nodal_stimulation_current = 40.                           # not used
    # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
    variables.vm_value_stimulated = 40.

elif "Aliev_Panfilov_Razumova_Titin" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
    # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
    variables.mappings = {
        ("parameter", 0): ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
        ("parameter", 1): ("constant", "Razumova/l_hs"),        # parameter 1 is constant 11 = fiber stretch λ
        # parameter 2 is constant 12 = fiber contraction velocity \dot{λ}
        ("parameter", 2): ("constant", "Razumova/rel_velo"),
        ("connectorSlot", "vm"): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
        # expose algebraic 4 = γ to the operator splitting
        ("connectorSlot", "stress"): ("algebraic", "Razumova/ActiveStress"),
        # expose algebraic 5 = α to the operator splitting
        ("connectorSlot", "alpha"): ("algebraic", "Razumova/Activation"),
    }
    # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
    variables.parameters_initial_values = [0, 1, 0]
    variables.nodal_stimulation_current = 40.                           # not used
    # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
    variables.vm_value_stimulated = 40.


# callback functions
# --------------------------
def get_motor_unit_no(fiber_no):
    return int(variables.fiber_distribution[fiber_no % len(variables.fiber_distribution)] - 1)


def get_diffusion_prefactor(fiber_no, mu_no):
    diffusion_prefactor = variables.get_conductivity(
        fiber_no, mu_no) / (variables.get_am(fiber_no, mu_no) * variables.get_cm(fiber_no, mu_no))
    # print("diffusion_prefactor: {}/({}*{}) = {}".format(variables.get_conductivity(fiber_no, mu_no), variables.get_am(fiber_no, mu_no), variables.get_cm(fiber_no, mu_no), diffusion_prefactor))
    return diffusion_prefactor


def fiber_gets_stimulated(fiber_no, frequency, current_time):
    """
    determine if fiber fiber_no gets stimulated at simulation time current_time
    """

    # determine motor unit
    alpha = 1.0   # 0.8
    mu_no = (int)(get_motor_unit_no(fiber_no) * alpha)

    # determine if fiber fires now
    index = int(np.round(current_time * frequency))
    n_firing_times = np.size(variables.firing_times, 0)

    # if variables.firing_times[index % n_firing_times, mu_no] == 1:
    # print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), variables.firing_times[index % n_firing_times, mu_no], "true" if variables.firing_times[index % n_firing_times, mu_no] == 1 else "false"))

    return variables.firing_times[index % n_firing_times, mu_no] == 1


def set_parameters(n_nodes_global, time_step_no, current_time, parameters, dof_nos_global, fiber_no):

    # determine if fiber gets stimulated at the current time
    is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, variables.stimulation_frequency, current_time)

    # determine nodes to stimulate (center node, left and right neighbour)
    innervation_zone_width_n_nodes = variables.innervation_zone_width * 100  # 100 nodes per cm
    # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    innervation_node_global = int(n_nodes_global / 2)
    nodes_to_stimulate_global = [innervation_node_global]
    if innervation_node_global > 0:
        nodes_to_stimulate_global.insert(0, innervation_node_global - 1)
    if innervation_node_global < n_nodes_global - 1:
        nodes_to_stimulate_global.append(innervation_node_global + 1)

    # stimulation value
    if is_fiber_gets_stimulated:
        stimulation_current = variables.nodal_stimulation_current
    else:
        stimulation_current = 0.

    first_dof_global = dof_nos_global[0]
    last_dof_global = dof_nos_global[-1]

    for node_no_global in nodes_to_stimulate_global:
        if first_dof_global <= node_no_global <= last_dof_global:
            # get local no for global no (1D)
            dof_no_local = node_no_global - first_dof_global
            parameters[dof_no_local] = stimulation_current

# callback function that can set parameters, i.e. stimulation current


def set_specific_parameters(n_nodes_global, time_step_no, current_time, parameters, fiber_no):

    # determine if fiber gets stimulated at the current time
    is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, variables.stimulation_frequency, current_time)

    # determine nodes to stimulate (center node, left and right neighbour)
    innervation_zone_width_n_nodes = variables.innervation_zone_width * 100  # 100 nodes per cm
    # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    innervation_node_global = int(n_nodes_global / 2)
    nodes_to_stimulate_global = [innervation_node_global]

    for k in range(10):
        if innervation_node_global - k >= 0:
            nodes_to_stimulate_global.insert(0, innervation_node_global - k)
        if innervation_node_global + k <= n_nodes_global - 1:
            nodes_to_stimulate_global.append(innervation_node_global + k)

    # stimulation value
    if is_fiber_gets_stimulated:
        stimulation_current = 40.
    else:
        stimulation_current = 0.

    for node_no_global in nodes_to_stimulate_global:
        parameters[(node_no_global, 0)] = stimulation_current   # key: ((x,y,z),nodal_dof_index)

# callback function that can set states, i.e. prescribed values for stimulation


def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):

    # print("call set_specific_states at time {}".format(current_time))

    # determine if fiber gets stimulated at the current time
    is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, variables.stimulation_frequency, current_time)

    if is_fiber_gets_stimulated:
        # determine nodes to stimulate (center node, left and right neighbour)
        innervation_zone_width_n_nodes = variables.innervation_zone_width * 100  # 100 nodes per cm
        # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
        innervation_node_global = int(n_nodes_global / 2)
        nodes_to_stimulate_global = [innervation_node_global]
        if innervation_node_global > 0:
            nodes_to_stimulate_global.insert(0, innervation_node_global - 1)
        if innervation_node_global < n_nodes_global - 1:
            nodes_to_stimulate_global.append(innervation_node_global + 1)
        # if rank_no == 0:
        #  print("t: {}, stimulate fiber {} at nodes {}".format(current_time, fiber_no, nodes_to_stimulate_global))

        for node_no_global in nodes_to_stimulate_global:
            states[(node_no_global, 0, 0)] = variables.vm_value_stimulated   # key: ((x,y,z),nodal_dof_index,state_no)

# callback function for artifical stress values, instead of monodomain


def set_stress_values(n_dofs_global, n_nodes_global_per_coordinate_direction,
                      time_step_no, current_time, values, global_natural_dofs, fiber_no):
    # n_dofs_global:       (int) global number of dofs in the mesh where to set the values
    # n_nodes_global_per_coordinate_direction (list of ints)   [mx, my, mz] number of global nodes in each coordinate direction.
    #                       For composite meshes, the values are only for the first submesh, for other meshes sum(...) equals n_dofs_global
    # time_step_no:        (int)   current time step number
    # current_time:        (float) the current simulation time
    # values:              (list of floats) all current local values of the field variable, if there are multiple components, they are stored in struct-of-array memory layout
    #                       i.e. [point0_component0, point0_component1, ... pointN_component0, point1_component0, point1_component1, ...]
    #                       After the call, these values will be assigned to the field variable.
    # global_natural_dofs  (list of ints) for every local dof no. the dof no. in global natural ordering
    # additional_argument: The value of the option "additionalArgument", can be any Python object.

    # loop over nodes in fiber
    for local_dof_no in range(len(values)):
        # get the global no. of the current dof
        global_dof_no = global_natural_dofs[local_dof_no]

        n_nodes_per_fiber = n_nodes_global_per_coordinate_direction[0]

        k = global_dof_no
        N = n_nodes_per_fiber

        if k > N / 2:
            k = N / 2 - k
        else:
            k = k - N / 2

        values[local_dof_no] = 0.1 * np.sin((current_time / 100 + 0.2 * k / N + 0.1 *
                                            fiber_no / variables.n_fibers_total) * 2 * np.pi) ** 2


# load MU distribution and firing times
variables.fiber_distribution = np.genfromtxt(variables.fiber_distribution_file, delimiter=" ")
variables.firing_times = np.genfromtxt(variables.firing_times_file)

# for debugging output show when the first 20 fibers will fire
if rank_no == 0 and not variables.disable_firing_output:
    print("Debugging output about fiber firing: Taking input from file \"{}\"".format(variables.firing_times_file))
    import timeit
    t_start = timeit.default_timer()

    first_stimulation_info = []

    n_firing_times = np.size(variables.firing_times, 0)
    for fiber_no_index in range(variables.n_fibers_total):
        if fiber_no_index % 100 == 0:
            t_algebraic = timeit.default_timer()
            if t_algebraic - t_start > 100:
                print("Note: break after {}/{} fibers ({:.0f}%) because it already took {:.3f}s".format(fiber_no_index,
                      variables.n_fibers_total, 100.0 * fiber_no_index / (variables.n_fibers_total - 1.), t_algebraic - t_start))
                break

        first_stimulation = None
        for current_time in np.linspace(0, 1. / variables.stimulation_frequency * n_firing_times, n_firing_times):
            if fiber_gets_stimulated(fiber_no_index, variables.stimulation_frequency, current_time):
                first_stimulation = current_time
                break
        mu_no = get_motor_unit_no(fiber_no_index)
        first_stimulation_info.append([fiber_no_index, mu_no, first_stimulation])

    first_stimulation_info.sort(
        key=lambda x: 1e6 +
        1e-6 *
        x[1] +
        1e-12 *
        x[0] if x[2] is None else x[2] +
        1e-6 *
        x[1] +
        1e-12 *
        x[0])

    print("First stimulation times")
    print("    Time  MU fibers")
    n_stimulated_mus = 0
    n_not_stimulated_mus = 0
    stimulated_fibers = []
    last_time = 0
    last_mu_no = first_stimulation_info[0][1]
    for stimulation_info in first_stimulation_info:
        mu_no = stimulation_info[1]
        fiber_no = stimulation_info[0]
        if mu_no == last_mu_no:
            stimulated_fibers.append(fiber_no)
        else:
            if last_time is not None:
                if len(stimulated_fibers) > 10:
                    print("{:8.2f} {:3} {} (only showing first 10, {} total)".format(
                        last_time, last_mu_no, str(stimulated_fibers[0:10]), len(stimulated_fibers)))
                else:
                    print("{:8.2f} {:3} {}".format(last_time, last_mu_no, str(stimulated_fibers)))
                n_stimulated_mus += 1
            else:
                if len(stimulated_fibers) > 10:
                    print("  never stimulated: MU {:3}, fibers {} (only showing first 10, {} total)".format(
                        last_mu_no, str(stimulated_fibers[0:10]), len(stimulated_fibers)))
                else:
                    print("  never stimulated: MU {:3}, fibers {}".format(last_mu_no, str(stimulated_fibers)))
                n_not_stimulated_mus += 1
            stimulated_fibers = [fiber_no]

        last_time = stimulation_info[2]
        last_mu_no = mu_no

    print("stimulated MUs: {}, not stimulated MUs: {}".format(n_stimulated_mus, n_not_stimulated_mus))

    t_end = timeit.default_timer()
    print("duration of assembling this list: {:.3f} s\n".format(t_end - t_start))

# compute partitioning
if rank_no == 0:
    if n_ranks != variables.n_subdomains_x * variables.n_subdomains_y * variables.n_subdomains_z:
        print(
            "\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(
                n_ranks,
                variables.n_subdomains_x,
                variables.n_subdomains_y,
                variables.n_subdomains_z,
                variables.n_subdomains_x *
                variables.n_subdomains_y *
                variables.n_subdomains_z))
        quit()

# n_fibers_per_subdomain_* is already set

####################################
# set Dirichlet BC for the flow problem

n_points_3D_mesh_linear_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x)
                                       for subdomain_coordinate_x in range(variables.n_subdomains_x)])
n_points_3D_mesh_linear_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y)
                                       for subdomain_coordinate_y in range(variables.n_subdomains_y)])
n_points_3D_mesh_linear_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z)
                                       for subdomain_coordinate_z in range(variables.n_subdomains_z)])
n_points_3D_mesh_linear_global = n_points_3D_mesh_linear_global_x * \
    n_points_3D_mesh_linear_global_y * n_points_3D_mesh_linear_global_z

n_points_3D_mesh_quadratic_global_x = 2 * n_points_3D_mesh_linear_global_x - 1
n_points_3D_mesh_quadratic_global_y = 2 * n_points_3D_mesh_linear_global_y - 1
n_points_3D_mesh_quadratic_global_z = 2 * n_points_3D_mesh_linear_global_z - 1

# set boundary conditions for the elasticity
[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
[nx, ny, nz] = variables.meshes["3Dmesh_quadratic"]["nElements"]

variables.fiber_mesh_names = [mesh_name for mesh_name in variables.meshes.keys() if "MeshFiber" in mesh_name]

# set Dirichlet BC at top nodes for linear elasticity problem, fix muscle at top
variables.elasticity_dirichlet_bc = {}
if False:
    for j in range(my):
        for i in range(mx):
            variables.elasticity_dirichlet_bc[(mz - 1) * mx * my + j * mx + i] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# fix muscle at bottom
if False:
    k = 0
    for j in range(my):
        for i in range(mx):
            variables.elasticity_dirichlet_bc[k * mx * my + j * mx + i] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Neumann BC at top nodes, traction upwards
# k = nz-1
# variables.elasticity_neumann_bc = [{"element": k*nx*ny + j*nx + i, "constantVector": [0.0,0.0,10.0], "face": "2+"} for j in range(ny) for i in range(nx)]
variables.elasticity_neumann_bc = []

# with open("mesh","w") as f:
#  f.write(str(variables.meshes["3Dmesh_quadratic"]))
