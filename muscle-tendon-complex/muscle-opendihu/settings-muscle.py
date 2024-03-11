# Multiple 1D fibers (monodomain) with 3D dynamic mooney rivlin with active contraction term, on biceps geometry

import sys
import os
import timeit
import argparse
import importlib
import distutils.util

# set title of terminal
title = "muscle"
print('\33]0;{}\a'.format(title), end='', flush=True)

# parse rank arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path, 'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there
# file create_partitioned_meshes_for_settings with helper functions about own subdomain
from create_partitioned_meshes_for_settings import *

# if first argument contains "*.py", it is a custom variable definition file, load these values
if ".py" in sys.argv[0]:
    variables_path_and_filename = sys.argv[0]
    variables_path, variables_filename = os.path.split(variables_path_and_filename)  # get path and filename
    # add the directory of the variables file to python path
    sys.path.insert(0, os.path.join(script_path, variables_path))
    # remove the ".py" extension to get the name of the module
    variables_module, _ = os.path.splitext(variables_filename)

    if rank_no == 0:
        print("Loading variables from \"{}\".".format(variables_path_and_filename))

    custom_variables = importlib.import_module(variables_module,
                                               package=variables_filename)    # import variables module
    variables.__dict__.update(custom_variables.__dict__)
    sys.argv = sys.argv[1:]     # remove first argument, which now has already been parsed
else:
    if rank_no == 0:
        print("Error: no variables file was specified, e.g:\n ./biceps_contraction ../settings_biceps_contraction.py ramp.py")
    exit(0)

# -------------- begin user parameters ----------------
variables.output_timestep_3D = 50  # [ms] output timestep of mechanics
variables.output_timestep_fibers = 50  # [ms] output timestep of fibers
# -------------- end user parameters ----------------

# define command line arguments


def mbool(x): return bool(distutils.util.strtobool(x))   # function to parse bool arguments


parser = argparse.ArgumentParser(description='muscle')
parser.add_argument(
    '--scenario_name',
    help='The name to identify this run in the log.',
    default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3, help='Number of subdomains in x,y,z direction.', type=int)
parser.add_argument(
    '--n_subdomains_x',
    '-x',
    help='Number of subdomains in x direction.',
    type=int,
    default=variables.n_subdomains_x)
parser.add_argument(
    '--n_subdomains_y',
    '-y',
    help='Number of subdomains in y direction.',
    type=int,
    default=variables.n_subdomains_y)
parser.add_argument(
    '--n_subdomains_z',
    '-z',
    help='Number of subdomains in z direction.',
    type=int,
    default=variables.n_subdomains_z)
parser.add_argument(
    '--diffusion_solver_type',
    help='The solver for the diffusion.',
    default=variables.diffusion_solver_type,
    choices=[
        "gmres",
        "cg",
        "lu",
        "gamg",
        "richardson",
        "chebyshev",
        "cholesky",
        "jacobi",
        "sor",
         "preonly"])
parser.add_argument(
    '--diffusion_preconditioner_type',
    help='The preconditioner for the diffusion.',
    default=variables.diffusion_preconditioner_type,
    choices=[
        "jacobi",
        "sor",
        "lu",
        "ilu",
        "gamg",
         "none"])
parser.add_argument(
    '--potential_flow_solver_type',
    help='The solver for the potential flow (non-spd matrix).',
    default=variables.potential_flow_solver_type,
    choices=[
        "gmres",
        "cg",
        "lu",
        "gamg",
        "richardson",
        "chebyshev",
        "cholesky",
        "jacobi",
        "sor",
         "preonly"])
parser.add_argument(
    '--potential_flow_preconditioner_type',
    help='The preconditioner for the potential flow.',
    default=variables.potential_flow_preconditioner_type,
    choices=[
        "jacobi",
        "sor",
        "lu",
        "ilu",
        "gamg",
         "none"])
parser.add_argument(
    '--paraview_output',
    help='Enable the paraview output writer.',
    default=variables.paraview_output,
    action='store_true')
parser.add_argument(
    '--adios_output',
    help='Enable the MegaMol/ADIOS output writer.',
    default=variables.adios_output,
    action='store_true')
parser.add_argument(
    '--fiber_file',
    help='The filename of the file that contains the fiber data.',
    default=variables.fiber_file)
parser.add_argument(
    '--fiber_distribution_file',
    help='The filename of the file that contains the MU firing times.',
    default=variables.fiber_distribution_file)
parser.add_argument(
    '--firing_times_file',
    help='The filename of the file that contains the cellml model.',
    default=variables.firing_times_file)
parser.add_argument(
    '--end_time',
    '--tend',
    '-t',
    help='The end simulation time.',
    type=float,
    default=variables.end_time)
parser.add_argument(
    '--output_timestep',
    help='The timestep for writing outputs.',
    type=float,
    default=variables.output_timestep)
parser.add_argument('--dt_0D', help='The timestep for the 0D model.', type=float, default=variables.dt_0D)
parser.add_argument('--dt_1D', help='The timestep for the 1D model.', type=float, default=variables.dt_1D)
parser.add_argument(
    '--dt_splitting',
    help='The timestep for the splitting.',
    type=float,
    default=variables.dt_splitting)
parser.add_argument(
    '--dt_3D',
    help='The timestep for the 3D model, i.e. dynamic solid mechanics.',
    type=float,
    default=variables.dt_3D)
parser.add_argument(
    '--disable_firing_output',
    help='Disables the initial list of fiber firings.',
    default=variables.disable_firing_output,
    action='store_true')
parser.add_argument(
    '--enable_coupling',
    help='Enables the precice coupling.',
    type=mbool,
    default=variables.enable_coupling)
parser.add_argument('--v', help='Enable full verbosity in c++ code')
parser.add_argument('-v', help='Enable verbosity level in c++ code', action="store_true")
parser.add_argument('-vmodule', help='Enable verbosity level for given file in c++ code')
parser.add_argument('-pause', help='Stop at parallel debugging barrier', action="store_true")

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)


variables.n_subdomains = variables.n_subdomains_x * variables.n_subdomains_y * variables.n_subdomains_z

# automatically initialize partitioning if it has not been set
if n_ranks != variables.n_subdomains:

    # create all possible partitionings to the given number of ranks
    optimal_value = n_ranks**(1 / 3)
    possible_partitionings = []
    for i in range(1, n_ranks + 1):
        for j in range(1, n_ranks + 1):
            if i * j <= n_ranks and n_ranks % (i * j) == 0:
                k = (int)(n_ranks / (i * j))
                performance = (k - optimal_value)**2 + (j - optimal_value)**2 + 1.1 * (i - optimal_value)**2
                possible_partitionings.append([i, j, k, performance])

    # if no possible partitioning was found
    if len(possible_partitionings) == 0:
        if rank_no == 0:
            print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {} and no automatic partitioning could be done.\n\n\033[0m".format(
                n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x * variables.n_subdomains_y * variables.n_subdomains_z))
        quit()

    # select the partitioning with the lowest value of performance which is the best
    lowest_performance = possible_partitionings[0][3] + 1
    for i in range(len(possible_partitionings)):
        if possible_partitionings[i][3] < lowest_performance:
            lowest_performance = possible_partitionings[i][3]
            variables.n_subdomains_x = possible_partitionings[i][0]
            variables.n_subdomains_y = possible_partitionings[i][1]
            variables.n_subdomains_z = possible_partitionings[i][2]

# output information of run
if rank_no == 0:
    print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name,
          variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
    print("dt_0D:           {:0.0e}, diffusion_solver_type:      {}".format(
        variables.dt_0D, variables.diffusion_solver_type))
    print("dt_1D:           {:0.0e}, potential_flow_solver_type: {}".format(
        variables.dt_1D, variables.potential_flow_solver_type))
    print("dt_splitting:    {:0.0e}, emg_solver_type:            {}, emg_initial_guess_nonzero: {}".format(
        variables.dt_splitting, variables.emg_solver_type, variables.emg_initial_guess_nonzero))
    print("dt_3D:           {:0.0e}, paraview_output: {}".format(variables.dt_3D, variables.paraview_output))
    print("output_timestep: {:0.0e}  stimulation_frequency: {} 1/ms = {} Hz".format(variables.output_timestep,
          variables.stimulation_frequency, variables.stimulation_frequency * 1e3))
    print("fiber_file:              {}".format(variables.fiber_file))
    print("cellml_file:             {}".format(variables.cellml_file))
    print("fiber_distribution_file: {}".format(variables.fiber_distribution_file))
    print("firing_times_file:       {}".format(variables.firing_times_file))
    print("********************************************************************************")

    print("prefactor: sigma_eff/(Am*Cm) = {} = {} / ({}*{})".format(variables.Conductivity /
          (variables.Am * variables.Cm), variables.Conductivity, variables.Am, variables.Cm))

    # start timer to measure duration of parsing of this script
    t_start_script = timeit.default_timer()

# initialize all helper variables
from helper import *

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

if False:
    for subdomain_coordinate_y in range(variables.n_subdomains_y):
        for subdomain_coordinate_x in range(variables.n_subdomains_x):

            print("subdomain (x{},y{}) ranks: {} n fibers in subdomain: x{},y{}".format(subdomain_coordinate_x, subdomain_coordinate_y,
                                                                                        list(
                                                                                            range(
                                                                                                subdomain_coordinate_y *
                                                                                                variables.n_subdomains_x +
                                                                                                subdomain_coordinate_x,
                                                                                                n_ranks,
                                                                                                variables.n_subdomains_x *
                                                                                                variables.n_subdomains_y)),
                                                                                        n_fibers_in_subdomain_x(subdomain_coordinate_x), n_fibers_in_subdomain_y(subdomain_coordinate_y)))

            for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)):
                for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)):
                    print("({},{}) n instances: {}".format(fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y,
                                                           n_fibers_in_subdomain_x(subdomain_coordinate_x) * n_fibers_in_subdomain_y(subdomain_coordinate_y)))


# define the config dict
config = {
    "scenarioName": variables.scenario_name,
    "logFormat": "csv",
    # output file of a diagram that shows data connection between solvers
    "solverStructureDiagramFile": "out/muscle_solver_structure.txt",
    "mappingsBetweenMeshesLogFile": "out/muscle_mappings_between_meshes.txt",  # log file of when mappings between meshes occur
    "Meshes": variables.meshes,
    "MappingsBetweenMeshes": variables.mappings_between_meshes,
    "Solvers": {
        "diffusionTermSolver": {  # solver for the implicit timestepping scheme of the diffusion time step
            "maxIterations": 1e4,
            "relativeTolerance": 1e-10,
            "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual
            "solverType": variables.diffusion_solver_type,
            "preconditionerType": variables.diffusion_preconditioner_type,
            "dumpFilename": "",   # "out/dump_"
            "dumpFormat": "matlab",
        },
        "potentialFlowSolver": {  # solver for the initial potential flow, that is needed to estimate fiber directions for the bidomain equation
            "relativeTolerance": 1e-10,
            "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual
            "maxIterations": 1e4,
            "solverType": variables.potential_flow_solver_type,
            "preconditionerType": variables.potential_flow_preconditioner_type,
            "dumpFilename": "",
            "dumpFormat": "matlab",
        },
        "mechanicsSolver": {   # solver for the dynamic mechanics problem
            "relativeTolerance": 1e-10,           # 1e-10 relative tolerance of the linear solver
            "absoluteTolerance": 1e-10,          # 1e-10 absolute tolerance of the residual of the linear solver
            "solverType": "preonly",      # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
            "preconditionerType": "lu",           # type of the preconditioner
            "maxIterations": 1e4,           # maximum number of iterations in the linear solver
            "snesMaxFunctionEvaluations": 1e8,    # maximum number of function iterations
            "snesMaxIterations": 140,            # maximum number of iterations in the nonlinear solver
            "snesRelativeTolerance": 1e-5,        # relative tolerance of the nonlinear solver
            "snesAbsoluteTolerance": 1e-5,        # absolute tolerance of the nonlinear solver
            "snesLineSearchType": "l2",           # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
            "snesRebuildJacobianFrequency": 3,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again
            "dumpFilename": "",            # dump system matrix and right hand side after every solve
            "dumpFormat": "matlab",      # default, ascii, matlab
        }
    },
    "PreciceAdapter": {        # precice adapter for muscle
        "timeStepOutputInterval": 100,                        # interval in which to display current timestep and time in console
        "timestepWidth": 1,                          # coupling time step width, must match the value in the precice config
        # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
        "couplingEnabled": variables.enable_coupling,
        "preciceConfigFilename": variables.precice_config_file,    # the preCICE configuration file
        # name of the own precice participant, has to match the name given in the precice xml config file
        "preciceParticipantName": "Muscle",
        "scalingFactor": 1,                          # a factor to scale the exchanged data, prior to communication
        # if the output writers should be called only after a time window of
        # precice is complete, this means the timestep has converged
        "outputOnlyConvergedTimeSteps": True,
        "preciceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
            {
                "preciceMeshName": "Muscle-Bottom-Mesh",         # precice name of the 2D coupling mesh
                "face": "2-",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
            },
            {
                "preciceMeshName": "Muscle-Top-A-Mesh",           # precice name of the 2D coupling mesh
                "face": "2+",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
            },
            {
                "preciceMeshName": "Muscle-Top-B-Mesh",           # precice name of the 2D coupling mesh
                "face": "2+",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
            },
        ],
        "preciceData": [
            {
                # mode is one of "read-displacements-velocities", "read-traction",
                # "write-displacements-velocities", "write-traction"
                "mode": "read-displacements-velocities",
                # name of the precice coupling surface mesh, as given in the precice xml settings file
                "preciceMeshName": "Muscle-Bottom-Mesh",
                # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
                "displacementsName": "Displacement",
                # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
                "velocitiesName": "Velocity",
            },
            {
                # mode is one of "read-displacements-velocities", "read-traction",
                # "write-displacements-velocities", "write-traction"
                "mode": "read-displacements-velocities",
                # name of the precice coupling surface mesh, as given in the precice xml settings file
                "preciceMeshName": "Muscle-Top-A-Mesh",
                # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
                "displacementsName": "Displacement",
                # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
                "velocitiesName": "Velocity",
            },
            {
                # mode is one of "read-displacements-velocities", "read-traction",
                # "write-displacements-velocities", "write-traction"
                "mode": "read-displacements-velocities",
                # name of the precice coupling surface mesh, as given in the precice xml settings file
                "preciceMeshName": "Muscle-Top-B-Mesh",
                # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
                "displacementsName": "Displacement",
                # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
                "velocitiesName": "Velocity",
            },
            {
                # mode is one of "read-displacements-velocities", "read-traction",
                # "write-displacements-velocities", "write-traction"
                "mode": "write-traction",
                # name of the precice coupling surface mesh, as given in the precice xml settings
                "preciceMeshName": "Muscle-Bottom-Mesh",
                # name of the traction "data", i.e. field variable, as given in the precice xml settings file
                "tractionName": "Traction",
            },
            {
                # mode is one of "read-displacements-velocities", "read-traction",
                # "write-displacements-velocities", "write-traction"
                "mode": "write-traction",
                # name of the precice coupling surface mesh, as given in the precice xml settings
                "preciceMeshName": "Muscle-Top-A-Mesh",
                # name of the traction "data", i.e. field variable, as given in the precice xml settings file
                "tractionName": "Traction",
            },
            {
                # mode is one of "read-displacements-velocities", "read-traction",
                # "write-displacements-velocities", "write-traction"
                "mode": "write-traction",
                # name of the precice coupling surface mesh, as given in the precice xml settings
                "preciceMeshName": "Muscle-Top-B-Mesh",
                # name of the traction "data", i.e. field variable, as given in the precice xml settings file
                "tractionName": "Traction",
            },
        ],

        "Coupling": {
            "timeStepWidth": variables.dt_3D,  # 1e-1
            "logTimeStepWidthAsKey": "dt_3D",
            "durationLogKey": "duration_total",
            "timeStepOutputInterval": 1,
            "endTime": variables.end_time,
            # transfer gamma to MuscleContractionSolver, the receiving slots are λ, λdot, γ
            "connectedSlotsTerm1To2": {1: 2},
            "connectedSlotsTerm2To1": None,       # transfer nothing back
            "Term1": {        # monodomain, fibers
                "MultipleInstances": {
                    "logKey": "duration_subdomains_xy",
                    "ranksAllComputedInstances": list(range(n_ranks)),
                    "nInstances": variables.n_subdomains_xy,
                    "instances":
                    [{
                        "ranks": list(range(subdomain_coordinate_y * variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x * variables.n_subdomains_y)),

                        # this is for the actual model with fibers
                        "StrangSplitting": {
                            # "numberTimeSteps": 1,
                            "timeStepWidth": variables.dt_splitting,  # 1e-1
                            "logTimeStepWidthAsKey": "dt_splitting",
                            "durationLogKey": "duration_monodomain",
                            "timeStepOutputInterval": 100,
                            "endTime": variables.dt_splitting,
                            # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
                            "connectedSlotsTerm1To2": [0, 1, 2],
                            "connectedSlotsTerm2To1": [0, None, 2],   # transfer the same back, this avoids data copy

                            "Term1": {      # CellML, i.e. reaction term of Monodomain equation
                                "MultipleInstances": {
                                    "logKey": "duration_subdomains_z",
                                    "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x) * n_fibers_in_subdomain_y(subdomain_coordinate_y),
                                    "instances":
                                    [{
                                        # these rank nos are local nos to the outer instance of MultipleInstances,
                                        # i.e. from 0 to number of ranks in z direction
                                        "ranks": list(range(variables.n_subdomains_z)),
                                        "Heun": {
                                            "timeStepWidth": variables.dt_0D,                         # timestep width of 0D problem
                                            # key under which the time step width will be written to the log file
                                            "logTimeStepWidthAsKey": "dt_0D",
                                            "durationLogKey": "duration_0D",                           # log key of duration for this solver
                                            "timeStepOutputInterval": 1e4,                                     # how often to print the current timestep
                                            "initialValues": [],                                      # no initial values are specified
                                            "dirichletBoundaryConditions": {},                                      # no Dirichlet boundary conditions are specified
                                            # filename for a vtp file that contains the Dirichlet boundary condition
                                            # nodes and their values, set to None to disable
                                            "dirichletOutputFilename": None,
                                            # the boundary conditions and initial values would be given as global
                                            # numbers
                                            "inputMeshIsGlobal": True,
                                            "checkForNanInf": True,                                    # abort execution if the solution contains nan or inf values
                                            "nAdditionalFieldVariables": 0,                                       # number of additional field variables
                                            "additionalSlotNames": "",

                                            "CellML": {
                                                "modelFilename": variables.cellml_file,                          # input C++ source file or cellml XML file
                                                # "statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                                                "statesInitialValues": variables.states_initial_values,                # initial values for new_slow_TK
                                                # if the equilibrium values of the states should be computed before the
                                                # simulation starts
                                                "initializeStatesToEquilibrium": False,
                                                # if initializeStatesToEquilibrium is enable, the timestep width to use to
                                                # solve the equilibrium equation
                                                "initializeStatesToEquilibriumTimestepWidth": 1e-4,

                                                # optimization parameters
                                                # "vc", "simd", "openmp" type of generated optimizated source file
                                                "optimizationType": "vc",
                                                # if optimizationType is "vc", whether the exponential function exp(x)
                                                # should be approximate by (1+x/n)^n with n=1024
                                                "approximateExponentialFunction": True,
                                                # compiler flags used to compile the optimized model code
                                                "compilerFlags": "-fPIC -O3 -march=native -Wno-deprecated-declarations -shared ",
                                                # if optimizationType is "openmp", the maximum number of threads to use.
                                                # Default value 0 means no restriction.
                                                "maximumNumberOfThreads": 0,

                                                # stimulation callbacks
                                                # "libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                                                # "setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                                                # "setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                                                # callback function that sets states like Vm, activation can be
                                                # implemented by using this method and directly setting Vm values, or by
                                                # using setParameters/setSpecificParameters
                                                "setSpecificStatesFunction": set_specific_states,
                                                # "setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                                                "setSpecificStatesCallInterval": 0,                                                               # 0 means disabled
                                                # set_specific_states should be called variables.stimulation_frequency
                                                # times per ms
                                                "setSpecificStatesCallFrequency": variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),
                                                # random value to add or substract to setSpecificStatesCallFrequency every
                                                # stimulation, this is to add random jitter to the frequency
                                                "setSpecificStatesFrequencyJitter": variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no),
                                                # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                                                "setSpecificStatesRepeatAfterFirstCall": 0.01,
                                                # [ms] first time when to call setSpecificStates
                                                "setSpecificStatesCallEnableBegin": variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),
                                                # last argument that will be passed to the callback functions
                                                # set_specific_states, set_specific_parameters, etc.
                                                "additionalArgument": fiber_no,

                                                # parameters to the cellml model
                                                # mappings between parameters and algebraics/constants and between
                                                # outputConnectorSlots and states, algebraics or parameters, they are
                                                # defined in helper.py
                                                "mappings": variables.mappings,
                                                # [0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                                                "parametersInitialValues": variables.parameters_initial_values,

                                                # reference to the fiber mesh
                                                "meshName": "MeshFiber_{}".format(fiber_no),
                                                # a file that will contain the times of stimulations
                                                "stimulationLogFilename": "out/stimulation.log",
                                            },
                                            "OutputWriter": [
                                                {"format": "Paraview",
                                                 "outputInterval": 1,
                                                 "filename": "out/" + variables.scenario_name + "/0D_states({},{})".format(fiber_in_subdomain_coordinate_x,
                                                                                                                           fiber_in_subdomain_coordinate_y),
                                                 "binary": True,
                                                 "fixedFormat": False,
                                                 "combineFiles": True,
                                                 "fileNumbering": "incremental"}
                                            ] if variables.states_output else []

                                        },
                                    } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                                        for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                                        for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                                        for motor_unit_no in [get_motor_unit_no(fiber_no)]],
                                }
                            },
                            "Term2": {     # Diffusion
                                "MultipleInstances": {
                                    "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x) * n_fibers_in_subdomain_y(subdomain_coordinate_y),
                                    "instances":
                                    [{
                                        # these rank nos are local nos to the outer instance of MultipleInstances,
                                        # i.e. from 0 to number of ranks in z direction
                                        "ranks": list(range(variables.n_subdomains_z)),
                                        "CrankNicolson": {
                                            "initialValues": [],                                      # no initial values are given
                                            # "numberTimeSteps":            1,
                                            "timeStepWidth": variables.dt_1D,                         # timestep width for the diffusion problem
                                            "timeStepWidthRelativeTolerance": 1e-10,
                                            # key under which the time step width will be written to the log file
                                            "logTimeStepWidthAsKey": "dt_1D",
                                            "durationLogKey": "duration_1D",                           # log key of duration for this solver
                                            "timeStepOutputInterval": 1e4,                                     # how often to print the current timestep
                                            # old Dirichlet BC that are not used in FastMonodomainSolver: {0:
                                            # -75.0036, -1: -75.0036},
                                            "dirichletBoundaryConditions": {},
                                            # filename for a vtp file that contains the Dirichlet boundary condition
                                            # nodes and their values, set to None to disable
                                            "dirichletOutputFilename": None,
                                            "inputMeshIsGlobal": True,                                    # initial values would be given as global numbers
                                            "solverName": "diffusionTermSolver",                   # reference to the linear solver
                                            # number of additional field variables that will be written to the output
                                            # file, here for stress
                                            "nAdditionalFieldVariables": 2,
                                            "additionalSlotNames": ["stress", "activation"],
                                            "checkForNanInf": True,                                    # abort execution if the solution contains nan or inf values

                                            "FiniteElementMethod": {
                                                "inputMeshIsGlobal": True,
                                                "meshName": "MeshFiber_{}".format(fiber_no),
                                                "solverName": "diffusionTermSolver",
                                                # resolves to Conductivity / (Am * Cm)
                                                "prefactor": get_diffusion_prefactor(fiber_no, motor_unit_no),
                                                "slotName": None,
                                            },
                                            "OutputWriter": [
                                                # {"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
                                                # {"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                                                # {"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                                                # {"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
                                            ]
                                        },
                                    } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                                        for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                                        for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                                        for motor_unit_no in [get_motor_unit_no(fiber_no)]],
                                    "OutputWriter": [
                                        {"format": "Paraview",
                                         "outputInterval": int(1 / variables.dt_splitting * variables.output_timestep_fibers),
                                         "filename": "out/fibers",
                                         "binary": True,
                                         "fixedFormat": False,
                                         "combineFiles": True,
                                         "fileNumbering": "incremental"}
                                    ],
                                },
                            },
                        },

                        # this is for biceps_contraction_no_cell, i.e. PrescribedValues instead of fibers
                        "GodunovSplitting": {   # this splitting scheme is only needed to replicate the solver structure as with the fibers
                            "timeStepWidth": variables.dt_3D,
                            "logTimeStepWidthAsKey": "dt_splitting",
                            "durationLogKey": "duration_prescribed_values",
                            "timeStepOutputInterval": 100,
                            "endTime": variables.dt_3D,
                            "connectedSlotsTerm1To2": [],
                            "connectedSlotsTerm2To1": [],   # transfer the same back, this avoids data copy

                            "Term1": {
                                "MultipleInstances": {
                                    "logKey": "duration_subdomains_z",
                                    "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x) * n_fibers_in_subdomain_y(subdomain_coordinate_y),
                                    "instances":
                                    [{
                                        # these rank nos are local nos to the outer instance of MultipleInstances,
                                        # i.e. from 0 to number of ranks in z direction
                                        "ranks": list(range(variables.n_subdomains_z)),
                                        "PrescribedValues": {
                                            # reference to the fiber mesh
                                            "meshName": "MeshFiber_{}".format(fiber_no),
                                            "numberTimeSteps": 1,             # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set the reaction term only once per time step
                                            "timeStepOutputInterval": 20,            # if the time step should be written to console, a value > 10 produces no output
                                            "slotNames": [],            # names of the data connector slots

                                            # a list of field variables that will get values assigned in every
                                            # timestep, by the provided callback function
                                            "fieldVariables1": [
                                                {"name": "Vm", "callback": None},
                                                {"name": "stress", "callback": set_stress_values},
                                            ],
                                            "fieldVariables2": [],
                                            # a custom argument to the fieldVariables callback functions, this will be
                                            # passed on as the last argument
                                            "additionalArgument": fiber_no,

                                            "OutputWriter": [
                                                {"format": "Paraview",
                                                 "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_fibers),
                                                 "filename": "out/prescribed_fibers",
                                                 "binary": True,
                                                 "fixedFormat": False,
                                                 "combineFiles": True,
                                                 "fileNumbering": "incremental"}
                                            ]
                                        },
                                    } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                                        for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                                        for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                                        for motor_unit_no in [get_motor_unit_no(fiber_no)]],

                                    # "OutputWriter" : variables.output_writer_fibers,
                                    "OutputWriter": [
                                        {"format": "Paraview",
                                         "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_fibers),
                                         "filename": "out/fibers",
                                         "binary": True,
                                         "fixedFormat": False,
                                         "combineFiles": True,
                                         "fileNumbering": "incremental"}
                                    ]
                                }
                            },

                            # term2 is unused, it is needed to be similar to the actual fiber solver structure
                            "Term2": {}
                        }

                    } if (subdomain_coordinate_x, subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y) else None
                        for subdomain_coordinate_y in range(variables.n_subdomains_y)
                        for subdomain_coordinate_x in range(variables.n_subdomains_x)]
                },
                # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
                "fiberDistributionFile": variables.fiber_distribution_file,
                "firingTimesFile": variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
                # only compute fibers after they have been stimulated for the first time
                "onlyComputeIfHasBeenStimulated": True,
                # optimization where states that are close to their equilibrium will not be computed again
                "disableComputationWhenStatesAreCloseToEquilibrium": True,
                "valueForStimulatedPoint": variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set
                # range where the neuromuscular junction is located around the center,
                # relative to fiber length. The actual position is draws randomly from the
                # interval [0.5-s/2, 0.5+s/2) with s being this option. 0 means sharply at
                # the center, 0.1 means located approximately at the center, but it can
                # vary 10% in total between all fibers.
                "neuromuscularJunctionRelativeSize": 0.1,
            },
            "Term2": {        # solid mechanics
                "MuscleContractionSolver": {
                    "numberTimeSteps": 1,                         # only use 1 timestep per interval
                    "timeStepOutputInterval": 100,                       # do not output time steps
                    "Pmax": variables.pmax,            # maximum PK2 active stress
                    # if the factor f_l(λ_f) modeling the force-length relation (as in
                    # Heidlauf2013) should be multiplied. Set to false if this relation is
                    # already considered in the CellML model.
                    "enableForceLengthRelation": variables.enable_force_length_relation,
                    # scaling factor for the output of the lambda dot slot, i.e. the
                    # contraction velocity. Use this to scale the unit-less quantity to, e.g.,
                    # micrometers per millisecond for the subcellular model.
                    "lambdaDotScalingFactor": variables.lambda_dot_scaling_factor,
                    "slotNames": ["lambda", "ldot", "gamma", "T"],   # names of the data connector slots
                    "OutputWriter": [
                        {"format": "Paraview",
                         "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D),
                         "filename": "out/muscle_3D",
                         "binary": True,
                         "fixedFormat": False,
                         "onlyNodalValues": True,
                         "combineFiles": True,
                         "fileNumbering": "incremental"},
                    ],
                    "mapGeometryToMeshes": [],                        # the mesh names of the meshes that will get the geometry transferred
                    "dynamic": True,                      # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem

                    # the actual solid mechanics solver, this is either
                    # "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on
                    # the value of "dynamic"
                    "DynamicHyperelasticitySolver": {
                        "timeStepWidth": variables.dt_3D,           # time step width
                        "durationLogKey": "nonlinear",               # key to find duration of this solver in the log file
                        "timeStepOutputInterval": 1,                         # how often the current time step should be printed to console

                        "materialParameters": variables.material_parameters,  # material parameters of the Mooney-Rivlin material
                        "density": variables.rho,             # density of the material
                        # scaling factor for displacements, only set to sth. other than 1 only to
                        # increase visual appearance for very small displacements
                        "displacementsScalingFactor": 1.0,
                        # log file where residual norm values of the nonlinear solver will be written
                        "residualNormLogFilename": "muscle_log_residual_norm.txt",
                        # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
                        "useAnalyticJacobian": True,
                        # whether to use the numerically computed jacobian matrix in the nonlinear
                        # solver (slow), only works with non-nested matrices, if both numeric and
                        # analytic are enable, it uses the analytic for the preconditioner and the
                        # numeric as normal jacobian
                        "useNumericJacobian": False,

                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
                        "dumpDenseMatlabVariables": False,
                        # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables
                        # all all three true, the analytic and numeric jacobian matrices will get
                        # compared to see if there are programming errors for the analytic
                        # jacobian

                        # mesh
                        "inputMeshIsGlobal": True,                     # the mesh is given locally
                        "meshName": "3Dmesh_quadratic",        # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
                        # fiber meshes that will be used to determine the fiber direction, for
                        # multidomain there are no fibers so this would be empty list
                        "fiberMeshNames": variables.fiber_mesh_names,
                        # "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system

                        # solving
                        # name of the nonlinear solver configuration, it is defined under
                        # "Solvers" at the beginning of this config
                        "solverName": "mechanicsSolver",
                        # "loadFactors":                [0.25, 0.66, 1.0],                # load factors for every timestep
                        "loadFactors": [],                        # no load factors, solve problem directly
                        # a threshold for the load factor, when to abort the solve of the current
                        # time step. The load factors are adjusted automatically if the nonlinear
                        # solver diverged. If the progression between two subsequent load factors
                        # gets smaller than this value, the solution is aborted.
                        "loadFactorGiveUpThreshold": 4e-2,
                        "nNonlinearSolveCalls": 1,                         # how often the nonlinear solve should be repeated

                        # boundary and initial conditions
                        # the initial Dirichlet boundary conditions that define values for
                        # displacements u and velocity v
                        "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,
                        # Neumann boundary conditions that define traction forces on surfaces of elements
                        "neumannBoundaryConditions": variables.elasticity_neumann_bc,
                        # if the given Neumann boundary condition values under
                        # "neumannBoundaryConditions" are total forces instead of surface loads
                        # and therefore should be scaled by the surface area of all elements where
                        # Neumann BC are applied
                        "divideNeumannBoundaryConditionValuesByTotalArea": True,
                        # function that updates the dirichlet BCs while the simulation is running
                        "updateDirichletBoundaryConditionsFunction": None,
                        # every which step the update function should be called, 1 means every time step
                        "updateDirichletBoundaryConditionsFunctionCallInterval": 1,

                        # the initial values for the displacements, vector of values for every
                        # node [[node1-x,y,z], [node2-x,y,z], ...]
                        "initialValuesDisplacements": [[0.0, 0.0, 0.0] for _ in range(mx * my * mz)],
                        # the initial values for the velocities, vector of values for every node
                        # [[node1-x,y,z], [node2-x,y,z], ...]
                        "initialValuesVelocities": [[0.0, 0.0, 0.0] for _ in range(mx * my * mz)],
                        # if the initial values for the dynamic nonlinear problem should be
                        # computed by extrapolating the previous displacements and velocities
                        "extrapolateInitialGuess": True,
                        "constantBodyForce": variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity

                        # filename for a vtp file that contains the Dirichlet boundary condition
                        # nodes and their values, set to None to disable
                        "dirichletOutputFilename": "out/muscle_dirichlet_boundary_conditions",
                        # filename of a log file that will contain the total (bearing) forces and
                        # moments at the top and bottom of the volume
                        "totalForceLogFilename": "out/muscle_force.csv",
                        "totalForceLogOutputInterval": 10,                                  # output interval when to write the totalForceLog file
                        # global element nos of the bottom elements used to compute the total
                        # forces in the log file totalForceLogFilename
                        "totalForceBottomElementNosGlobal": [j * nx + i for j in range(ny) for i in range(nx)],
                        # global element nos of the top elements used to compute the total forces
                        # in the log file totalForceTopElementsGlobal
                        "totalForceTopElementNosGlobal": [(nz - 1) * ny * nx + j * nx + i for j in range(ny) for i in range(nx)],


                        # define which file formats should be written
                        # 1. main output writer that writes output files using the quadratic
                        # elements function space. Writes displacements, velocities and PK2
                        # stresses.
                        "OutputWriter": [

                            # Paraview files
                            # {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},

                            # Python callback function "postprocess"
                            # {"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
                        ],
                        # 2. additional output writer that writes also the hydrostatic pressure
                        "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
                            "OutputWriter": [
                                # {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                            ]
                        },
                        # 3. additional output writer that writes virtual work terms
                        "dynamic": {    # output of the dynamic solver, has additional virtual work values
                            "OutputWriter": [   # output files for displacements function space (quadratic elements)
                                # {"format": "Paraview", "outputInterval": int(output_interval/dt), "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                                {"format": "Paraview",
                                 "outputInterval": int(1. / variables.dt_3D * variables.output_timestep_3D),
                                 "filename": "out/muscle_virtual_work",
                                 "binary": True,
                                 "fixedFormat": False,
                                 "onlyNodalValues": True,
                                 "combineFiles": True,
                                 "fileNumbering": "incremental"},
                            ],
                        },
                        # 4. output writer for debugging, outputs files after each load increment,
                        # the geometry is not changed but u and v are written
                        "LoadIncrements": {
                            "OutputWriter": [
                                # {"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                            ]
                        },
                    }
                }
            }
        }
    }
}


# stop timer and calculate how long parsing lasted
if rank_no == 0:
    t_stop_script = timeit.default_timer()
    print("Python config parsed in {:.1f}s.".format(t_stop_script - t_start_script))
