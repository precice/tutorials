
using Downloads: download
using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

function initial_condition_channel(x, t, equations::CompressibleEulerEquations2D)
    rho = 1.0
    rho_v1 = 0.2
    rho_v2 = 0.0
    rho_e = 10.0
    return SVector(rho, rho_v1, rho_v2, rho_e)
end

initial_condition = initial_condition_channel

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

###############################################################################
# Create mesh

coordinates_min = (0.0, 0.0) # minimum coordinates (min(x), min(y))
coordinates_max = (6.0, 2.0) # maximum coordinates (max(x), max(y))

# Create a uniformly refined mesh
trees_per_dimension = (6, 2)
mesh = P4estMesh(trees_per_dimension,
                 polydeg=1, initial_refinement_level=1,
                 coordinates_min=coordinates_min, coordinates_max=coordinates_max,
                 periodicity=(true, true))
                 # periodicity=(false, false))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    # boundary_conditions=Dict(
                                    #   :all => BoundaryConditionDirichlet(initial_condition)
                                    # )
                                   )


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

# save_solution = SaveSolutionCallback(interval=100,
#                                      save_initial_solution=true,
#                                      save_final_solution=true,
#                                      solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        # save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary
