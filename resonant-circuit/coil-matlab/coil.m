clear; close all; clc;

% Initialize and configure preCICE
interface = precice.Participant("Coil", "../precice-config.xml", 0, 1);

% Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
meshName ="Coil-Mesh";
vertex_ID = interface.setMeshVertex(meshName, [0 0]);

% Data IDs
dataNameI = "Current";
dataNameU = "Voltage";

% Simulation parameters and initial condition
L = 1;                      % Inductance
t0 = 0;                     % Initial simulation time
t_max = 10;                 % End simulation time
Io = 1;                     % Initial voltage
phi = 0;                    % Phase of the signal

I0 = Io*cos(phi);           % Initial condition for I

f_I = @(t, I, U) U/L;       % Time derivative of I

if interface.requiresInitialData()
    interface.writeData(meshName, DataNameI, vertex_ID, I0);
end
interface.initialize();
dt = interface.getMaxTimeStepSize();

% Initialize simulation
U0 = interface.readData(meshName, dataNameU, vertex_ID, 0); % Initial voltage
t_vec = t0;                 % Vector of time


% Start simulation
t = t0 + dt;
while interface.isCouplingOngoing()

    % Record checkpoint if necessary
    if interface.requiresWritingCheckpoint()
        I0_checkpoint = I0;
        U0_checkpoint = U0;
    end

    % Make Simulation Step
    [t_ode, I_ode] = ode45(@(t, y) f_I(t, y, U0), [t0 t], I0);
    I0 = I_ode(end);

    % Exchange data
    interface.writeData(meshName, dataNameI, vertex_ID, I0);
    interface.advance(dt);

    % Recover checkpoint if not converged, else finish time step
    if interface.requiresReadingCheckpoint()
        I0 = I0_checkpoint;
        U0 = U0_checkpoint;
    else
        dt = interface.getMaxTimeStepSize();
        U0 = interface.readData(meshName, dataNameU, vertex_ID, dt);
        t0 = t;
        t = t0 + dt;
    end

end

% Stop coupling
interface.finalize();
