clear; close all; clc;

% Initialize and configure preCICE
interface = precice.Participant("ParticipantU", "../precice-config.xml", 0, 1);

% Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
meshName ="MeshU";
vertex_ID = interface.setMeshVertex(meshName, [0 0]);

% Data IDs
dataNameI = "I";
dataNameU = "U";

% Simulation parameters and initial condition
C = 2;                      % Capacitance
L = 1;                      % Inductance
t0 = 0;                     % Initial simulation time
t_max = 10;                 % End simulation time
Io = 1;                     % Initial current
phi = 0;                    % Phase of the signal

w0 = 1/sqrt(L*C);           % Resonant frequency
I0 = Io*cos(phi);           % Initial condition for I
U0 = -w0*L*Io*sin(phi);     % Initial condition for U

f_U = @(t, U, I) -I/C;      % Time derivative of U

% Initialize simulation
if interface.requiresInitialData()
    interface.writeData(meshName, DataNameU, vertex_ID, U0);
end
interface.initialize();
dt = interface.getMaxTimeStepSize();

% Start simulation
t = t0 + dt;
while interface.isCouplingOngoing()

    % Record checkpoint if necessary
    if interface.requiresWritingCheckpoint()
        I0_checkpoint = I0;
        U0_checkpoint = U0;
    end

    % Make Simulation Step
    [t_ode, U_ode] = ode45(@(t, y) f_U(t, y, I0), [t0 t], U0);
    U0 = U_ode(end);

    % Exchange data
    interface.writeData(meshName, dataNameU, vertex_ID, U0);
    interface.advance(dt);

    % Recover checkpoint if not converged, else finish time step
    if interface.requiresReadingCheckpoint()
        I0 = I0_checkpoint;
        U0 = U0_checkpoint;
    else
        dt = interface.getMaxTimeStepSize();
        I0 = interface.readData(meshName, dataNameI, vertex_ID, dt);
        t0 = t;
        t = t0 + dt;
    end

end

% Stop coupling
interface.finalize();
