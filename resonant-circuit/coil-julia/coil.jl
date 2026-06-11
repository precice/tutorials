using PreCICE
using DifferentialEquations

# Initialize and configure preCICE
participant = PreCICE.createParticipant("Coil", "../precice-config.xml", 0, 1)

# Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
mesh_name = "Coil-Mesh"

dimensions = PreCICE.getMeshDimensions(mesh_name)

vertex_ids = PreCICE.setMeshVertices(mesh_name, zeros((1,dimensions)))

let
    # Data IDs
    read_data_name = "Voltage"
    write_data_name = "Current"

    # Simulation parameters and initial condition
    C = 2                      # Capacitance of the capacitor
    L = 1                      # Inductance of the coil
    t0 = 0                     # Initial simulation time
    Io = 1                     # Initial current
    phi = 0                    # Phase of the signal

    w0 = 1 / sqrt(L * C)       # Resonant frequency
    I0 = Io * cos(phi)         # Initial condition for I

    # to estimate cost
    f_evals = 0

    function f_I(u, p, t)
        f_evals += 1
        (dt, L, mesh_name, read_data_name, vertex_ids) = p
        U = PreCICE.readData(mesh_name, read_data_name, vertex_ids, dt)
        return U[1] / L  # Time derivative of I; ODE determining capacitor
    end

    # Initialize simulation
    if PreCICE.requiresInitialData()
        PreCICE.writeData(mesh_name, write_data_name, vertex_ids, I0)
    end
    PreCICE.initialize()

    solver_dt = PreCICE.getMaxTimeStepSize()

    # Start simulation
    t = t0
    I0_checkpoint = I0
    t_checkpoint = t
    while PreCICE.isCouplingOngoing()

        # Record checkpoint if necessary
        if PreCICE.requiresWritingCheckpoint()
            I0_checkpoint = I0
            t_checkpoint = t
        end

        # Make Simulation Step
        precice_dt = PreCICE.getMaxTimeStepSize()
        dt = min(precice_dt, solver_dt)
        t_span = (t, t + dt)
        params = (dt, L, mesh_name, read_data_name, vertex_ids)
        prob = ODEProblem(f_I, I0, t_span, params)
        # https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)

        # Exchange data
        evals = max(length(sol.t), 3)  # at least do 3 substeps to allow cubic interpolation
        for i in range(1,evals)
            I0 = sol(t_span[1] + i * dt / evals)
            PreCICE.writeData(mesh_name, write_data_name, vertex_ids, fill(I0, (1,1)))
            PreCICE.advance(dt / evals)
        end
        t = t + dt

        # Recover checkpoint if not converged
        if PreCICE.requiresReadingCheckpoint()
            I0 = I0_checkpoint
            t = t_checkpoint
        end
    end
    # Stop coupling
    PreCICE.finalize()
end