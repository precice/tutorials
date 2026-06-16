using PreCICE
using DifferentialEquations

# Initialize and configure preCICE
participant = PreCICE.createParticipant("Capacitor", "../precice-config.xml", 0, 1)

# Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
mesh_name = "Capacitor-Mesh"

dimensions = PreCICE.getMeshDimensions(mesh_name)

vertex_ids = PreCICE.setMeshVertices(mesh_name, zeros((1,dimensions)))

let
    # Data IDs
    read_data_name = "Current"
    write_data_name = "Voltage"

    # Simulation parameters and initial condition
    C = 2                         # Capacitance of the capacitor
    L = 1                         # Inductance of the coil
    t0 = 0                        # Initial simulation time
    t_max = 10                    # End simulation time
    Io = 1                        # Initial current
    phi = 0                       # Phase of the signal

    w0 = 1 / sqrt(L * C)          # Resonant frequency
    U0 = -w0 * L * Io * sin(phi)  # Initial condition for U

    function f_U(u, p, t)
        (dt, C, mesh_name, read_data_name, vertex_ids) = p
        I = PreCICE.readData(mesh_name, read_data_name, vertex_ids, dt)
        return -I[1] / C      # Time derivative of U
    end

    # Initialize simulation
    if PreCICE.requiresInitialData()
        PreCICE.writeData(mesh_name, write_data_name, vertex_ids, I0)
    end
    PreCICE.initialize()

    solver_dt = PreCICE.getMaxTimeStepSize()

    # Start simulation
    t = t0
    U0_checkpoint = U0
    t_checkpoint = t
    while PreCICE.isCouplingOngoing()

        # Record checkpoint if necessary
        if PreCICE.requiresWritingCheckpoint()
            U0_checkpoint = U0
            t_checkpoint = t
        end

        # Make Simulation Step
        precice_dt = PreCICE.getMaxTimeStepSize()
        dt = min(precice_dt, solver_dt)
        t_span = (t, t + dt)
        params = (dt, C, mesh_name, read_data_name, vertex_ids)
        prob = ODEProblem(f_U, U0, t_span, params)
        # https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)

        # Exchange data
        evals = max(length(sol.t), 3)  # at least do 3 substeps to allow cubic interpolation
        for i in range(1,evals)
            U0 = sol(t_span[1] + i * dt / evals)
            PreCICE.writeData(mesh_name, write_data_name, vertex_ids, fill(U0, (1,1)))
            PreCICE.advance(dt / evals)
        end
        t = t + dt

        # Recover checkpoint if not converged
        if PreCICE.requiresReadingCheckpoint()
            U0 = U0_checkpoint
            t = t_checkpoint
        end
    end
    # Stop coupling
    PreCICE.finalize()
end