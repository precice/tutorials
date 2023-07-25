import sys
import os
import itertools
import importlib

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])
var_file = sys.argv[0] if ".py" in sys.argv[0] else "variables.py"

# add folders to python path
script_path = os.path.dirname(os.path.abspath(__file__))
var_path = os.path.join(script_path, "../variables")
sys.path.insert(0, var_path)

# load variables file
scenario_name, _ = os.path.splitext(var_file)
variables = importlib.import_module(scenario_name, package=var_file)

# define config
config = {
  "scenarioName":                   scenario_name,

  "logFormat":                      "csv",
  "mappingsBetweenMeshesLogFile":   "out/" + scenario_name + "/mappings_between_meshes.txt",
  "solverStructureDiagramFile":     "out/" + scenario_name + "/solver_structure.txt",
  
  "Meshes":                         variables.meshes,
  "MappingsBetweenMeshes":          {},

  "Solvers": {
    "mechanicsSolver": {
      "solverType":                     "preonly",
      "preconditionerType":             "lu",
      "relativeTolerance":              1e-10,
      "absoluteTolerance":              1e-10,
      "maxIterations":                  1e4,
      "snesLineSearchType":             "l2",
      "snesRelativeTolerance":          1e-5,
      "snesAbsoluteTolerance":          1e-5,
      "snesMaxIterations":              10,
      "snesMaxFunctionEvaluations":     1e8,
      "snesRebuildJacobianFrequency":   5,
      "dumpFilename":                   "",
      "dumpFormat":                     "matlab"
    }
  },

  "RepeatedCall": {
    "timeStepWidth":            variables.dt_3D,
    "logTimeStepWidthAsKey":    "dt_3D",
    "durationLogKey":           "duration_3D",
    "endTime":                  variables.end_time,
    "timeStepOutputInterval":   1,

    "MuscleContractionSolver": {
      "Pmax":                         variables.pmax,
      "slotNames":                    ["lambda", "ldot", "gamma", "T"],
      "dynamic":                      True,
  
      "numberTimeSteps":              1,
      "timeStepOutputInterval":       100,
      "lambdaDotScalingFactor":       1,
      "enableForceLengthRelation":    True,
      "mapGeometryToMeshes":          [],
  
      "OutputWriter": [
        {
          "format":             "Paraview",
          "outputInterval":     int(1.0 / variables.dt_3D * variables.output_interval),
          "filename":           "out/" + scenario_name + "/muscle",
          "fileNumbering":      "incremental",
          "binary":             True,
          "fixedFormat":        False,
          "onlyNodalValues":    True,
          "combineFiles":       True
        }
      ],
  
      "DynamicHyperelasticitySolver": {
        "durationLogKey":         "duration_3D",
        "logTimeStepWidthAsKey":  "dt_3D",
        "numberTimeSteps":        1,
        "materialParameters":     variables.material_parameters,
        "density":                variables.rho,
        "timeStepOutputInterval": 1,
  
        "meshName":                   "mesh3D",
        "fiberDirectionInElement":    variables.fiber_direction,
        "inputMeshIsGlobal":          True,
        "fiberMeshNames":             [],
        "fiberDirection":             None,
  
        "solverName":                 "mechanicsSolver",
        "displacementsScalingFactor":  1.0,
        "useAnalyticJacobian":        True,
        "useNumericJacobian":         False,
        "dumpDenseMatlabVariables":   False,
        "loadFactorGiveUpThreshold":  1,
        "loadFactors":                [],
        "scaleInitialGuess":          False,
        "extrapolateInitialGuess":    True,
        "nNonlinearSolveCalls":       1,
  
        "dirichletBoundaryConditions":                            variables.dirichlet_bc,
        "neumannBoundaryConditions":                              variables.neumann_bc,
        "updateDirichletBoundaryConditionsFunction":              None,
        "updateDirichletBoundaryConditionsFunctionCallInterval":  1,
        "divideNeumannBoundaryConditionValuesByTotalArea":        True,
  
        "initialValuesDisplacements": [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
        "initialValuesVelocities":    [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
        "constantBodyForce":          (0, 0, 0),
  
        "dirichletOutputFilename":    "out/" + scenario_name + "/dirichlet_output",
        "residualNormLogFilename":    "out/" + scenario_name + "/residual_norm_log.txt",
        "totalForceLogFilename":      "out/" + scenario_name + "/total_force_log.txt",
  
        "OutputWriter": [],
        "pressure":       { "OutputWriter": [] },
        "dynamic":        { "OutputWriter": [] },
        "LoadIncrements": { "OutputWriter": [] }
      }
    }
  }
}
