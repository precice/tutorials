clear; close all; clc;

% Initialize and configure preCICE
interface = precice.Participant("ParticipantI", "../precice-config.xml", 0, 1);

% Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
meshName ="MeshI";
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

f_I = @(t, I, U) U/L;       % Time derivative of I

% Initialize simulation
I = I0;                     % Vector of I through time
U = U0;                     % Vector of U through time
t_vec = t0;                 % Vector of time
if interface.requiresInitialData()
    interface.writeData(meshName, DataNameI, vertex_ID, I0);
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
        U = [U U0];
        I = [I I0];
        t_vec = [t_vec, t];
        t0 = t;
        t = t0 + dt;
    end

end

% Stop coupling
interface.finalize();

% Analytical solution for comparison
I_an = Io*cos(w0*t_vec+phi);
U_an = -w0*L*Io*sin(w0*t_vec+phi);

% Make and save plot
figure(1)
subplot(2,1,1)
plot(t_vec, I_an)
hold on;
plot(t_vec, U_an)
ylim([-1,1])
legend('I', 'U')
title('Analytical')

subplot(2,1,2)
plot(t_vec, I)
hold on;
plot(t_vec, U)
ylim([-1,1])
title('Numerical')
legend('I', 'U')

save('outputs.mat', 'I', 'U')
saveas(gcf, 'Curves.png')
