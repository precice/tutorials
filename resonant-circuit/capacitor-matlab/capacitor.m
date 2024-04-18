clear; close all; clc;

% Initialize and configure preCICE
interface = precice.Participant("Capacitor", "../precice-config.xml", 0, 1);

% Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
meshName ="Capacitor-Mesh";
vertex_ID = interface.setMeshVertex(meshName, [0 0]);

% Data IDs
dataNameI = "Current";
dataNameU = "Voltage";

% Simulation parameters and initial condition
C = 2;                      % Capacitance
L = 1;                      % Inductance
t0 = 0;                     % Initial simulation time
t_max = 10;                 % End simulation time
Io = 1;                     % Initial voltage
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

I = I0;                     % Vector of I through time
U = U0;                     % Vector of U through time
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

% print maximal error over the whole time interval for the voltage and the
% current
error_I = max(abs(I_an-I));
error_U = max(abs(U_an-U));
fprintf('Error in I: %d \n', error_I);
fprintf('Error in U: %d \n', error_U);

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
