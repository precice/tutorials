# Transversely-isotropic Mooney Rivlin on a tendon geometry
# Note, this is not possible to be run in parallel because the fibers
# cannot be initialized without MultipleInstances class.
import sys
import os
import numpy as np
import sys

# set title of terminal
title = "tendon-top-a"
print('\33]0;{}\a'.format(title), end='', flush=True)

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path, 'variables'))

import variables
from create_partitioned_meshes_for_settings import *

# update material parameters
if (variables.tendon_material == "nonLinear"):
    c = 9.98                    # [N/cm^2=kPa]
    ca = 14.92                  # [-]
    ct = 14.7                   # [-]
    cat = 9.64                  # [-]
    ctt = 11.24                 # [-]
    mu = 3.76                   # [N/cm^2=kPa]
    k1 = 42.217e3               # [N/cm^2=kPa]
    k2 = 411.360e3              # [N/cm^2=kPa]

    variables.material_parameters = [c, ca, ct, cat, ctt, mu, k1, k2]

if (variables.tendon_material == "SaintVenantKirchoff"):
    # material parameters for Saint Venant-Kirchhoff material
    # https://www.researchgate.net/publication/230248067_Bulk_Modulus

    youngs_modulus = 7e4        # [N/cm^2 = 10kPa]
    shear_modulus = 3e4

    lambd = shear_modulus * (youngs_modulus - 2 * shear_modulus) / \
        (3 * shear_modulus - youngs_modulus)  # Lamé parameter lambda
    mu = shear_modulus       # Lamé parameter mu or G (shear modulus)

    variables.material_parameters = [lambd, mu]

# If the fiber geometry data should be loaded completely in the python
# script. If True, this reads the binary file and assigns the node
# positions in the config. If False, the C++ code will read the binary
# file and only extract the local node positions. This is more performant
# for highly parallel runs.
load_fiber_data = False

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

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
        sys.exit(-1)

# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 2

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z,
    variables.fiber_file, load_fiber_data,
    sampling_stride_x, sampling_stride_y, sampling_stride_z, True, True)

[variables.meshes,
    variables.own_subdomain_coordinate_x,
    variables.own_subdomain_coordinate_y,
    variables.own_subdomain_coordinate_z,
    variables.n_fibers_x,
    variables.n_fibers_y,
 variables.n_points_whole_fiber] = result

n_points_3D_mesh_linear_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x)
                                       for subdomain_coordinate_x in range(variables.n_subdomains_x)])
n_points_3D_mesh_linear_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y)
                                       for subdomain_coordinate_y in range(variables.n_subdomains_y)])
n_points_3D_mesh_linear_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z)
                                       for subdomain_coordinate_z in range(variables.n_subdomains_z)])
n_points_3D_mesh_linear_global = n_points_3D_mesh_linear_global_x * \
    n_points_3D_mesh_linear_global_y * n_points_3D_mesh_linear_global_z
nx = n_points_3D_mesh_linear_global_x - 1
ny = n_points_3D_mesh_linear_global_y - 1
nz = n_points_3D_mesh_linear_global_z - 1

node_positions = variables.meshes["3Dmesh_quadratic"]["nodePositions"]

# boundary conditions (for quadratic elements)
# --------------------------------------------
[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
[nx, ny, nz] = variables.meshes["3Dmesh_quadratic"]["nElements"]

# set Dirichlet BC, fix top end of tendon that is attached to the bone
variables.elasticity_dirichlet_bc = {}
k = mz - 1

# fix the whole x-y plane
for j in range(my):
    for i in range(mx):
        variables.elasticity_dirichlet_bc[k * mx * my + j * mx + i] = [0.0, 0.0, 0.0, None, None, None]

# set no Neumann BC
variables.elasticity_neumann_bc = []


config_hyperelasticity = {    # for both "HyperelasticitySolver" and "DynamicHyperelasticitySolver"
    "timeStepWidth": variables.dt_elasticity,      # time step width
    "endTime": variables.end_time,           # end time of the simulation time span
    "durationLogKey": "duration_mechanics",         # key to find duration of this solver in the log file
    "timeStepOutputInterval": 1,                            # how often the current time step should be printed to console

    "materialParameters": variables.material_parameters,  # material parameters of the Mooney-Rivlin material
    "density": variables.rho,                # density of the material
    # scaling factor for displacements, only set to sth. other than 1 only to
    # increase visual appearance for very small displacements
    "displacementsScalingFactor": 1.0,
    # log file where residual norm values of the nonlinear solver will be written
    "residualNormLogFilename": "out/tendon_top_a_log_residual_norm.txt",
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
    "meshName": "3Dmesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
    # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
    "inputMeshIsGlobal": True,

    "fiberMeshNames": [],                           # fiber meshes that will be used to determine the fiber direction
    # "fiberDirection":             [0,0,1],                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
    # if fiberMeshNames and fiberDirections are empty, directly set the
    # constant fiber direction, in element coordinate system
    "fiberDirectionInElement": [0, 0, 1],

    # nonlinear solver
    "relativeTolerance": 1e-10,                         # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance": 1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver
    "solverType": "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",                         # type of the preconditioner
    "maxIterations": 1e4,                          # maximum number of iterations in the linear solver
    "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
    "snesMaxIterations": 240,                           # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance": 1e-2,                         # relative tolerance of the nonlinear solver
    # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesLineSearchType": "l2",
    "snesAbsoluteTolerance": 1e-5,                         # absolute tolerance of the nonlinear solver
    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild,
    # 1 means rebuild every time the Jacobian is computed within a single
    # nonlinear solve, 2 means every second time the Jacobian is built etc. -2
    # means rebuild at next chance but then never again
    "snesRebuildJacobianFrequency": 5,

    # "dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
    "dumpFilename": "",                           # dump disabled
    "dumpFormat": "matlab",                     # default, ascii, matlab

    # "loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
    # "loadFactors":                [0.5, 1.0],                   # load factors for every timestep
    "loadFactors": [],                           # no load factors, solve problem directly
    # a threshold for the load factor, when to abort the solve of the current
    # time step. The load factors are adjusted automatically if the nonlinear
    # solver diverged. If the load factors get too small, it aborts the solve.
    "loadFactorGiveUpThreshold": 1e-3,
    "nNonlinearSolveCalls": 1,                            # how often the nonlinear solve should be called

    # boundary and initial conditions
    # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
    "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,
    # Neumann boundary conditions that define traction forces on surfaces of elements
    "neumannBoundaryConditions": variables.elasticity_neumann_bc,
    # if the given Neumann boundary condition values under
    # "neumannBoundaryConditions" are total forces instead of surface loads
    # and therefore should be scaled by the surface area of all elements where
    # Neumann BC are applied
    "divideNeumannBoundaryConditionValuesByTotalArea": False,
    # function that updates the dirichlet BCs while the simulation is running
    "updateDirichletBoundaryConditionsFunction": None,
    # every which step the update function should be called, 1 means every time step
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,

    # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "initialValuesDisplacements": [[0.0, 0.0, 0.0] for _ in range(mx * my * mz)],
    # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "initialValuesVelocities": [[0.0, 0.0, 0.0] for _ in range(mx * my * mz)],
    # if the initial values for the dynamic nonlinear problem should be
    # computed by extrapolating the previous displacements and velocities
    "extrapolateInitialGuess": False,
    "constantBodyForce": variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity

    # filename for a vtp file that contains the Dirichlet boundary condition
    # nodes and their values, set to None to disable
    "dirichletOutputFilename": "out/tendon_top_a_dirichlet_boundary_conditions_tendon_top_a",
    # filename of a log file that will contain the total (bearing) forces and
    # moments at the top and bottom of the volume
    "totalForceLogFilename": "out/tendon_top_a_force.csv",
    "totalForceLogOutputInterval": 10,                                  # output interval when to write the totalForceLog file
    # global element nos of the bottom elements used to compute the total forces in the log file totalForceLogFilename
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
        {"format": "Paraview",
         "outputInterval": int(1. / variables.dt_elasticity * variables.output_timestep_3D),
         "filename": "out/tendon_top_a",
         "binary": True,
         "fixedFormat": False,
         "onlyNodalValues": True,
         "combineFiles": True,
         "fileNumbering": "incremental"},

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
            {"format": "Paraview",
             "outputInterval": int(1. / variables.dt_elasticity * variables.output_timestep_3D),
             "filename": "out/tendon_top_a_virtual_work",
             "binary": True,
             "fixedFormat": False,
             "onlyNodalValues": True,
             "combineFiles": True,
             "fileNumbering": "incremental"},
            # {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
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

config = {
    "scenarioName": variables.scenario_name,      # scenario name to identify the simulation runs in the log file
    "logFormat": "csv",                        # "csv" or "json", format of the lines in the log file, csv gives smaller files
    # output file of a diagram that shows data connection between solvers
    "solverStructureDiagramFile": "out/tendon_top_a_solver_structure.txt",
    "mappingsBetweenMeshesLogFile": "out/tendon_top_a_mappings_between_meshes_log.txt",    # log file for mappings
    "Meshes": variables.meshes,

    "PreciceAdapter": {        # precice adapter for bottom tendon
        "timeStepOutputInterval": 100,                        # interval in which to display current timestep and time in console
        "timestepWidth": 1,                          # coupling time step width, must match the value in the precice config
        # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
        "couplingEnabled": True,
        "preciceConfigFilename": variables.precice_config_file,    # the preCICE configuration file
        # name of the own precice participant, has to match the name given in the precice xml config file
        "preciceParticipantName": "Tendon-Top-A",
        "scalingFactor": 1,                          # a factor to scale the exchanged data, prior to communication
        # if the output writers should be called only after a time window of
        # precice is complete, this means the timestep has converged
        "outputOnlyConvergedTimeSteps": True,
        "preciceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
            {
                "preciceMeshName": "Tendon-Top-A-Mesh",        # precice name of the 2D coupling mesh
                "face": "2-",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
            }
        ],
        "preciceData": [
            {
                # mode is one of "read-displacements-velocities", "read-traction",
                # "write-displacements-velocities", "write-traction"
                "mode": "write-displacements-velocities",
                # name of the precice coupling surface mesh, as given in the precice xml settings file
                "preciceMeshName": "Tendon-Top-A-Mesh",
                # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
                "displacementsName": "Displacement",
                # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
                "velocitiesName": "Velocity",
            },
            {
                # mode is one of "read-displacements-velocities", "read-traction",
                # "write-displacements-velocities", "write-traction"
                "mode": "read-traction",
                # name of the precice coupling surface mesh, as given in the precice xml settings
                "preciceMeshName": "Tendon-Top-A-Mesh",
                # name of the traction "data", i.e. field variable, as given in the precice xml settings file
                "tractionName": "Traction",
            }
        ],
        "HyperelasticitySolver": config_hyperelasticity,
        "DynamicHyperelasticitySolver": config_hyperelasticity,
    }
}
