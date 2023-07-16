import sys
import os
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
  "mappingsBetweenMeshesLogFile":   "out/" + scenario_name + "/fibers_mappings_between_meshes_log.txt",
  "solverStructureDiagramFile":     "out/" + scenario_name + "/fibers_solver_structure_diagram.txt",
  
  "Meshes":                         variables.meshes,
  "MappingsBetweenMeshes":          {},

  "Solvers": {
    "diffusionSolver": {
      "solverType":                     "cg",
      "preconditionerType":             "none",
      "relativeTolerance":              1e-10,
      "absoluteTolerance":              1e-10,
      "maxIterations":                  1e4,
      "dumpFilename":                   "",
      "dumpFormat":                     "matlab"
    }
  },

  "PreciceAdapterVolumeCoupling": {
    "preciceConfigFilename":        variables.precice_config,
    "preciceParticipantName":       "Fibers",
    "couplingEnabled":              True,
    "timeStepOutputInterval":       100,
    "timestepWidth":                variables.dt_3D,
    "scalingFactor":                1,
    "outputOnlyConvergedTimeSteps": True,

    "preciceData": [
      {
        "mode":             "read",
        "preciceDataName":  "Geometry",
        "preciceMeshName":  "FibersMesh",
        "opendihuMeshName": None,
        "slotName":         None,
        "isGeometryField":  True
      },
      {
        "mode":             "write",
        "preciceDataName":  "Gamma",
        "preciceMeshName":  "FibersMesh",
        "opendihuMeshName": None,
        "slotName":         "stress",
        "isGeometryField":  False
       }
    ],

    "MultipleInstances": {
      "ranksAllComputedInstances":  list(range(n_ranks)),
      "nInstances":                 1,

      "instances": [{
        "ranks": [0],

        "StrangSplitting": {
          "timeStepWidth":              variables.dt_splitting,
          "logTimeStepWidthAsKey":      "dt_splitting",
          "durationLogKey":             "duration_splitting",
          "timeStepOutputInterval":     100,
          "connectedSlotsTerm1To2":     None,
          "connectedSlotsTerm2To1":     None,

          "Term1": { # reaction term
            "MultipleInstances": {
              "nInstances": variables.fb_x * variables.fb_y,

              "instances": [{
                "ranks": [0],

                "Heun": {
                  "timeStepWidth":          variables.dt_0D,
                  "logTimeStepWidthAsKey":  "dt_0D",
                  "durationLogKey":         "duration_0D",
                  "timeStepOutputInterval": 100,

                  "initialValues":                  [],
                  "dirichletBoundaryConditions":    {},
                  "dirichletOutputFilename":        None,
                  "inputMeshIsGlobal":              True,
                  "checkForNanInf":                 False,
                  "nAdditionalFieldVariables":      0,
                  "additionalSlotNames":            [],
                  "OutputWriter":                   [],

                  "CellML": {
                    "modelFilename":            variables.input_dir + "hodgkin_huxley-razumova.cellml",
                    "meshName":                 "fiber{}".format(variables.get_fiber_no(fiber_x, fiber_y)),
                    "stimulationLogFilename":   "out/" + scenario_name + "/fibers_stimulation.log",

                    "statesInitialValues":                          [],
                    "initializeStatesToEquilibrium":                False,
                    "initializeStatesToEquilibriumTimeStepWidth":   1e-4,
                    "optimizationType":                             "vc",
                    "approximateExponentialFunction":               True,
                    "compilerFlags":                                "-fPIC -O3 -march=native -Wno-deprecated_declarations -shared",
                    "maximumNumberOfThreads":                       0,

                    "setSpecificStatesCallEnableBegin":         variables.specific_states_call_enable_begin,
                    "setSpecificStatesCallFrequency":           variables.specific_states_call_frequency,
                    "setSpecificStatesRepeatAfterFirstCall":    0.01,
                    "setSpecificStatesFrequencyJitter":         [0],
                    "setSpecificStatesCallInterval":            0,
                    "setSpecificStatesFunction":                None,
                    "additionalArgument":                       None,

                    "mappings": {
                      ("parameter", 0):                 "membrane/i_Stim",
                      ("parameter", 1):                 "Razumova/l_hs",
                      ("parameter", 2):                 ("constant", "Razumova/rel_velo"),
                      ("connectorSlot", "vm"):          "membrane/V",
                      ("connectorSlot", "stress"):      "Razumova/activestress",
                      ("connectorSlot", "alpha"):       "Razumova/activation",
                      ("connectorSlot", "lambda"):      "Razumova/l_hs",
                      ("connectorSlot", "ldot"):        "Razumova/rel_velo"
                    },
                    "parametersInitialValues": [0.0, 1.0, 0.0]
                  }
                }
              } for fiber_x in range(variables.fb_x) for fiber_y in range(variables.fb_y)]
            }
          },

          "Term2": { # diffusion term
            "MultipleInstances": {
              "nInstances": variables.fb_x * variables.fb_y,

              "OutputWriter": [
                {
                  "format":         "Paraview",
                  "outputInterval": int(1.0 / variables.dt_splitting * variables.output_interval),
                  "filename":       "out/" + scenario_name + "/fibers",
                  "fileNumbering":  "incremental",
                  "binary":         True,
                  "fixedFormat":    False,
                  "combineFiles":   True
                }
              ],
              
              "instances": [{
                "ranks": [0],

                "ImplicitEuler": {
                  "timeStepWidth":          variables.dt_1D,
                  "logTimeStepWidthAsKey":  "dt_1D",
                  "durationLogKey":         "duration_1D",
                  "timeStepOutputInterval": 100,

                  "nAdditionalFieldVariables":  4,
                  "additionalSlotNames":        ["stress", "alpha", "lambda", "ldot"],
                  
                  "solverName":                     "diffusionSolver",
                  "timeStepWidthRelativeTolerance": 1e-10,

                  "dirichletBoundaryConditions":    {},
                  "dirichletOutputFilename":        None,
                  "inputMeshIsGlobal":              True,
                  "checkForNanInf":                 False,
                  "OutputWriter":                   [],

                  "FiniteElementMethod": {
                    "meshName":             "fiber{}".format(variables.get_fiber_no(fiber_x, fiber_y)),
                    "inputMeshIsGlobal":    True,
                    "solverName":           "diffusionSolver",
                    "prefactor":            variables.diffusion_prefactor,
                    "slotName":             "vm"
                  }
                }
              } for fiber_x in range(variables.fb_x) for fiber_y in range(variables.fb_y)]
            }
          }
        }
      }]
    },

    "fiberDistributionFile":                                variables.fiber_distribution_file,
    "firingTimesFile":                                      variables.firing_times_file,
    "valueForStimulatedPoint":                              20.0,
    "onlyComputeIfHasBeenStimulated":                       True,
    "disableComputationWhenStatesAreCloseToEquilibrium":    True,
    "neuromuscularJunctionRelativeSize":                    0.1,
    "generateGPUSource":                                    True,
    "useSinglePrecision":                                   False
  }
}
