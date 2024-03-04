# Transversely-isotropic Mooney Rivlin on a tendon geometry
# Note, this is not possible to be run in parallel because the fibers cannot be initialized without MultipleInstances class.
import sys, os
import numpy as np
import pickle
import argparse
import sys
sys.path.insert(0, '.')
import variables              # file variables.py, defines default values for all parameters, you can set the parameters there
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain
import stl
from stl import mesh
import json

# set title of terminal
title = "tendon_bottom"
print('\33]0;{}\a'.format(title), end='', flush=True)

# material parameters
# --------------------
# quantities in mechanics unit system
variables.rho = 10          # [1e-4 kg/cm^3] 10 = density of the muscle (density of water)

# material parameters for Saint Venant-Kirchhoff material
# https://www.researchgate.net/publication/230248067_Bulk_Modulus

youngs_modulus = 7e4        # [N/cm^2 = 10kPa]  
shear_modulus = 3e4

#youngs_modulus*=1e-3
#shear_modulus*=1e-3

lambd = shear_modulus*(youngs_modulus - 2*shear_modulus) / (3*shear_modulus - youngs_modulus)  # Lamé parameter lambda
mu = shear_modulus       # Lamé parameter mu or G (shear modulus)

variables.material_parameters = [lambd, mu]

variables.constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
variables.force = 100.0           # [N] pulling force to the bottom 

variables.dt_elasticity = 1      # [ms] time step width for elasticity
variables.end_time      = 20000   # [ms] simulation time
variables.scenario_name = "tendon_bottom"
variables.is_bottom_tendon = True        # whether the tendon is at the bottom (negative z-direction), this is important for the boundary conditions
variables.output_timestep_3D = 50  # [ms] output timestep

# input mesh file
fiber_file = "../../../../input/left_biceps_brachii_tendon1.bin"        # bottom tendon
#fiber_file = "../../../../input/left_biceps_brachii_tendon2a.bin"
#fiber_file = "../../../../input/left_biceps_brachii_tendon2b.bin"
#fiber_file = "../../../../input/left_biceps_brachii_7x7fibers.bin"
#fiber_file = "../../../../input/left_biceps_brachii_7x7fibers.bin"

load_fiber_data = False             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# define command line arguments
parser = argparse.ArgumentParser(description='tendon')
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('-vmodule', help='ignore')

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

# partitioning
# ------------
# this has to match the total number of processes
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]

# compute partitioning
if rank_no == 0:
  if n_ranks != variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z:
    print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    sys.exit(-1)
    
# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 2

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, 
    fiber_file, load_fiber_data,
    sampling_stride_x, sampling_stride_y, sampling_stride_z, True, True)

[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result

n_points_3D_mesh_linear_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)])
n_points_3D_mesh_linear_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)])
n_points_3D_mesh_linear_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)])
n_points_3D_mesh_linear_global = n_points_3D_mesh_linear_global_x*n_points_3D_mesh_linear_global_y*n_points_3D_mesh_linear_global_z
nx = n_points_3D_mesh_linear_global_x-1
ny = n_points_3D_mesh_linear_global_y-1
nz = n_points_3D_mesh_linear_global_z-1

node_positions = variables.meshes["3Dmesh_quadratic"]["nodePositions"]

# boundary conditions (for quadratic elements)
# --------------------------------------------
[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
[nx, ny, nz] = variables.meshes["3Dmesh_quadratic"]["nElements"]

# set Dirichlet BC, fix x and y coordinates of the end of tendon that is attached to the bone
variables.elasticity_dirichlet_bc = {}
k = 0
  
# fix z value on the whole x-y-plane
for j in range(my):
  for i in range(mx):
    variables.elasticity_dirichlet_bc[k*mx*my + j*mx + i] = [0.0,0.0,None,None,None,None]
       
# set Neumann BC, set traction at the end of the tendon that is attached to the bone
k = 0
# start with 0 BC
variables.elasticity_neumann_bc = [{"element": k*nx*ny + j*nx + i, "constantVector": [0,0,0], "face": "2-"} for j in range(ny) for i in range(nx)]

def update_neumann_bc(t):

  # set new Neumann boundary conditions
  k = 0
  factor = min(1, t/100)   # for t ∈ [0,100] from 0 to 1
  elasticity_neumann_bc = [{
		"element": k*nx*ny + j*nx + i, 
		"constantVector": [0,0,-variables.force*factor], 		# force pointing to bottom
		"face": "2-",
    "isInReferenceConfiguration": True
  } for j in range(ny) for i in range(nx)]

  config = {
    "inputMeshIsGlobal": True,
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    "neumannBoundaryConditions": elasticity_neumann_bc,
  }
  print("prescribed pulling force to bottom: {}".format(variables.force*factor))
  return config

# update dirichlet boundary conditions to account for movement of humerus

current_ulna_force = 0
current_ulna_angle = 0

# global coordinates of tendon bottom
# global coordinates of elbow hinge
elbow_hinge_point = np.array([3.54436, 11.4571, -58.5607])
bottom_tendon_insertion_point = np.array([4.30, 14.81, -63.41])

vec = -elbow_hinge_point + bottom_tendon_insertion_point
angle_offset = np.arctan((bottom_tendon_insertion_point[2] - elbow_hinge_point[2]) / np.linalg.norm(vec))
rotation_axis = np.array([-1.5, 1, 0])
rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)   # normalize rotation axis
  
def rotation_matrix(angle):
  axis_x = rotation_axis[0]
  axis_y = rotation_axis[1]
  axis_z = rotation_axis[2]
  
  # compute rotation matrix
  rotation_matrix = np.array([
    [np.cos(angle) + axis_x**2*(1 - np.cos(angle)), 
     axis_x*axis_y*(1 - np.cos(angle)) - axis_z*np.sin(angle), 
     axis_x*axis_z*(1 - np.cos(angle)) + axis_y*np.sin(angle)],
    [axis_y*axis_x*(1 - np.cos(angle)) + axis_z*np.sin(angle),
     np.cos(angle) + axis_y**2*(1 - np.cos(angle)),
     axis_y*axis_z*(1 - np.cos(angle)) - axis_x*np.sin(angle)],
    [axis_z*axis_x*(1 - np.cos(angle)) - axis_y*np.sin(angle),
     axis_z*axis_y*(1 - np.cos(angle)) + axis_x*np.sin(angle),
     np.cos(angle) + axis_z**2*(1 - np.cos(angle))]])

  return rotation_matrix

call_count = 0
ulna_series_files = []
ulna_stl_filename = os.path.join(os.getcwd(),"cm_left_ulna.stl")
stl_mesh = None
try:
  stl_mesh = mesh.Mesh.from_file(ulna_stl_filename)
except: 
  print("Could not open file {} in directory {}".format(ulna_stl_filename, os.getcwd()))
  sys.exit(0)

# Function to update dirichlet boundary conditions over time, t.
# Only those entries can be updated that were also initially set.
def update_dirichlet_bc(t):
  global current_ulna_angle

  # determine parameter for the rotation of the tendon insertion point
  angle = current_ulna_angle
  rotation_point = elbow_hinge_point
  rotation_mat = rotation_matrix(angle)
  vertex = bottom_tendon_insertion_point
  
  # rotation vertex about rotation_axis by angle
  vertex = vertex - rotation_point
  vertex = rotation_mat.dot(vertex)
  vertex = vertex + rotation_point

  new_insertion_point = vertex
  offset = -bottom_tendon_insertion_point + new_insertion_point
  print("angle: {} deg, rot matrix: {}".format(angle*180/np.pi,rotation_matrix(angle)))
  print("old insertion point: {}, new insertion point: {}, offset: {}".format(bottom_tendon_insertion_point,new_insertion_point,offset))
  #offset[0] = 0
  #offset[1] = 0
  #offset[2] = 0

  # update dirichlet boundary conditions, set prescribed value to offset, do not constrain velocity
  for key in variables.elasticity_dirichlet_bc.keys():
    variables.elasticity_dirichlet_bc[key] = [offset[0],offset[1],offset[2],None,None,None]
  return variables.elasticity_dirichlet_bc
 
def store_rotated_ulna(t): 
  global current_ulna_angle, call_count, stl_mesh, ulna_series_files

  if np.isnan(current_ulna_angle):
    return

  # store rotated ulna
  call_count += 1
  output_interval = 50    # [ms] (because coupling timestep is 1ms)
  if call_count % output_interval == 0:
    out_triangles = []

    rotation_point = elbow_hinge_point
    rotation_mat = rotation_matrix(current_ulna_angle)

    for p in stl_mesh.points:
      # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

      # transform vertices
      vertex_list = []
      
      # apply rotation
      for vertex in [np.array(p[0:3]), np.array(p[3:6]), np.array(p[6:9])]:
        vertex = vertex - rotation_point
        vertex = rotation_mat.dot(vertex)
        vertex = vertex + rotation_point
        vertex_list.append(vertex)
      
      out_triangles += [vertex_list]

    # Create the mesh
    out_mesh = mesh.Mesh(np.zeros(len(out_triangles), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(out_triangles):
      out_mesh.vectors[i] = f
        
    out_mesh.update_normals()
    ulna_output_filename = "ulna_{:04d}.stl".format((int)(call_count/output_interval))
    ulna_output_path = os.path.join("out",ulna_output_filename)
    out_mesh.save(ulna_output_path)
    print("Saved file {}".format(ulna_output_path))
    ulna_series_files.append({"name": ulna_output_filename, "time": t})

    # save json series file
    ulna_series_filename = "out/ulna.stl.series"
    with open(ulna_series_filename, "w") as f:  
      data = {"file-series-version" : "1.0", "files" : ulna_series_files}
      json.dump(data, f, indent='\t')

forces = []
def callback_total_force(t, bearing_force_bottom, bearing_moment_bottom, bearing_force_top, bearing_moment_top):
  global current_ulna_angle
  # this callback functions gets the current total forces that are exerted by the muscle  

  current_ulna_force = bearing_force_bottom[2]

  # compute average of last 10 values
  forces.append(current_ulna_force)

  if len(forces) > 10:
    forces.pop(0)
  current_ulna_force = np.mean(forces)
  print("callback_total_force, t: {}, force: {}".format(t, current_ulna_force))
  
  # compute relation between force and angle of ulna
  if False:
    # positive angle = elbow flexion (ulna move downwards)
    min_ulna_angle = 0  # [deg]
    max_ulna_angle = -45   # [deg]
    force_factor = -current_ulna_force / 500
    force_factor = min(1.5, max(-0.5, force_factor))

    current_ulna_angle = min_ulna_angle + force_factor * (max_ulna_angle-min_ulna_angle)
    print("callback_total_force, t: {}, force: {}, factor: {}, angle: {}".format(t, current_ulna_force, force_factor,current_ulna_angle))
    current_ulna_angle *= np.pi/180   # convert from deg to rad



# Function to postprocess the output
# This function gets periodically called by the running simulation. 
# It provides all current variables for each node: geometry (position), u, v, stress, etc.
def postprocess(result):
  global current_ulna_angle

  result = result[0]
  # print result for debugging
  #print(result)
  
  # get current time
  current_time = result["currentTime"]
  timestep_no = result["timeStepNo"]

  # get number of nodes
  nx = result["nElementsLocal"][0]		# number of elements
  ny = result["nElementsLocal"][1]    # number of elements
  nz = result["nElementsLocal"][2]    # number of elements
  mx = 2*nx + 1  # number of nodes for quadratic elements
  my = 2*ny + 1
  mz = 2*nz + 1
  
  # parse variables
  field_variables = result["data"]
  
  #for f in field_variables:
  #  print(f["name"])
  
  # field_variables[0] is the geometry
  # field_variables[1] is the displacements u
  # etc., uncomment the above to see all field variables
  
  displacement_components = field_variables[1]["components"]
  
  # traction values contains the traction vector in reference configuration
  u1_values = displacement_components[0]["values"]   # displacement in x-direction
  u2_values = displacement_components[1]["values"]   # displacement in y-direction
  u3_values = displacement_components[2]["values"]   # displacement in z-direction
  u3_values_bottom = [u3_values[j*mx+i] for j in range(my) for i in range(mx)]
  z_displacement = np.mean(u3_values)

  #print("nx,ny: {},{}, mx,my: {},{}, {}={}".format(nx,ny,mx,my,mx*my*mz,len(u3_values)))

  # compute elbow angle from displacement of muscle in z direction
  vec = -elbow_hinge_point + bottom_tendon_insertion_point
  current_ulna_angle = -np.arcsin(z_displacement / np.linalg.norm(vec))
 
  print("z displacement: {} ({}), current_ulna_angle: {} deg".format(z_displacement, np.linalg.norm(vec), current_ulna_angle*180/np.pi))

  store_rotated_ulna(current_time)

config_hyperelasticity = {    # for both "HyperelasticitySolver" and "DynamicHyperelasticitySolver"
  "timeStepWidth":              variables.dt_elasticity,      # time step width 
  "endTime":                    variables.end_time,           # end time of the simulation time span    
  "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
  "timeStepOutputInterval":     1,                            # how often the current time step should be printed to console
  
  "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
  "density":                    variables.rho,                # density of the material
  "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
  "residualNormLogFilename":    "out/tendon_bottom_log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
  "useAnalyticJacobian":        True,                         # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
  "useNumericJacobian":         False,                        # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
    
  "dumpDenseMatlabVariables":   False,                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
  # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
  
  # mesh
  "preciceMeshName":                   "3Dmesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
  "inputMeshIsGlobal":          True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
  
  "fiberMeshNames":             [],                           # fiber meshes that will be used to determine the fiber direction
  #"fiberDirection":             [0,0,1],                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
  "fiberDirectionInElement":    [0,0,1],                      # if fiberMeshNames and fiberDirections are empty, directly set the constant fiber direction, in element coordinate system
      
  # nonlinear solver
  "relativeTolerance":          1e-10,                         # 1e-10 relative tolerance of the linear solver
  "absoluteTolerance":          1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver       
  "solverType":                 "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
  "preconditionerType":         "lu",                         # type of the preconditioner
  "maxIterations":              1e4,                          # maximum number of iterations in the linear solver
  "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
  "snesMaxIterations":          240,                           # maximum number of iterations in the nonlinear solver
  "snesRelativeTolerance":      1e-2,                         # relative tolerance of the nonlinear solver
  "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
  "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
  "snesRebuildJacobianFrequency": 5,                          # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
  
  #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
  "dumpFilename":               "",                           # dump disabled
  "dumpFormat":                 "matlab",                     # default, ascii, matlab
  
  #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
  #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
  "loadFactors":                [],                           # no load factors, solve problem directly
  "loadFactorGiveUpThreshold":  1e-3,                         # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
  "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
  
  # boundary and initial conditions
  "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
  "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
  "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
  "updateDirichletBoundaryConditionsFunction": None, #update_dirichlet_bc,   # function that updates the dirichlet BCs while the simulation is running
  "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # stide every which step the update function should be called, 1 means every time step
  "updateNeumannBoundaryConditionsFunction": update_neumann_bc,       # a callback function to periodically update the Neumann boundary conditions
  "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step 
 
  "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
  "constantBodyForce":           variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
  
  "dirichletOutputFilename":     "out/tendon_bottom_dirichlet_boundary_conditions",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
  "totalForceLogFilename":       "out/tendon_bottom_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
  "totalForceLogOutputInterval": 10,                                  # output interval when to write the totalForceLog file
  "totalForceBottomElementNosGlobal":  [j*nx + i for j in range(ny) for i in range(nx)],                  # global element nos of the bottom elements used to compute the total forces in the log file totalForceLogFilename
  "totalForceTopElementNosGlobal":     [(nz-1)*ny*nx + j*nx + i for j in range(ny) for i in range(nx)],   # global element nos of the top elements used to compute the total forces in the log file totalForceTopElementsGlobal
  "totalForceFunction":          None, #callback_total_force,                # callback function that gets the total force at bottom and top of the domain
  "totalForceFunctionCallInterval": 1,                                # how often the "totalForceFunction" is called
      
  # define which file formats should be written
  # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
  "OutputWriter" : [
    
    # Paraview files
    {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_3D), "filename": "out/tendon_bottom", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    
    # Python callback function "postprocess"
    {"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": "", "fileNumbering": "incremental"},
  ],
  # 2. additional output writer that writes also the hydrostatic pressure
  "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    ]
  },
  # 3. additional output writer that writes virtual work terms
  "dynamic": {    # output of the dynamic solver, has additional virtual work values 
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_3D), "filename": "out/tendon_bottom_virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    ],
  },
  # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
  "LoadIncrements": {   
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    ]
  },
}

config = {
  "scenarioName":                   variables.scenario_name,      # scenario name to identify the simulation runs in the log file
  "logFormat":                      "csv",                        # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "out/tendon_bottom_solver_structure.txt",       # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/tendon_bottom_mappings_between_meshes_log.txt",    # log file for mappings 
  "Meshes":                         variables.meshes,
  
  "PreciceAdapter": {        # precice adapter for bottom tendon
    "timeStepOutputInterval":   100,                        # interval in which to display current timestep and time in console
    "timestepWidth":            1,                          # coupling time step width, must match the value in the precice config
    "couplingEnabled":          True,                       # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
    "preciceConfigFilename":    "precice_config_muscle_dirichlet_tendon_neumann_implicit_coupling_multiple_tendons.xml",    # the preCICE configuration file
    "preciceParticipantName":   "TendonSolverBottom",       # name of the own precice participant, has to match the name given in the precice xml config file
    "scalingFactor":            1,                          # a factor to scale the exchanged data, prior to communication
    "outputOnlyConvergedTimeSteps": True,                   # if the output writers should be called only after a time window of precice is complete, this means the timestep has converged
    "preciceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
      {
        "preciceMeshName":      "Tendon-Bottom-Mesh",            # precice name of the 2D coupling mesh
        "face":                 "2+",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
      }
    ],
    "preciceData": [
      {
        "mode":                 "write-displacements-velocities",   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "Tendon-Bottom-Mesh",                    # name of the precice coupling surface mesh, as given in the precice xml settings file
        "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
        "velocitiesName":       "Velocity",                         # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
      },
      {
        "mode":                 "read-traction",                    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "Tendon-Bottom-Mesh",                    # name of the precice coupling surface mesh, as given in the precice xml settings 
        "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
      }
    ],
    "HyperelasticitySolver": config_hyperelasticity,
    "DynamicHyperelasticitySolver": config_hyperelasticity,
  }
}
