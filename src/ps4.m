close all; clear; clc

%% Problem 1(a)
% Recalculate moments of inertia
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

% Find Euler kinematics
tFinal = 60;
tStep = 0.01;
t = 0:tStep:tFinal;

eulerAngle0 = [0; 90; 0];
w0 = [0; 0; 1];
state0 = [eulerAngle0;w0];

tspan = 0:tStep:tFinal;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinEulerAngle(t,state,Ix,Iy,Iz), ...
    tspan,state0,options);

state = wrapTo360(rad2deg(state(:,:)));
w = state(:,4:6);

% Plot
figure()
plot(t,state(:,4:6),'LineWidth',1)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(1,'Images/ps4_problem1a_velocity.png')

figure()
plot(t,state(:,1:3),'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel(['Euler Angle [' char(176) ']'])
saveas(1,'Images/ps4_problem1a_angle.png')

%% Problem 1(b)
a = 7125.48662; % km
e = 0.0011650;
i = 98.40508; % degree
O = -19.61601; % degree
w_deg = 89.99764; % degree
nu = -89.99818; % degree

[~,y] = plotECI(a,e,i,O,w_deg,nu,t);
close all

w_RTN = zeros(size(w));

for n = 1:length(t)
    % Get rotation matrix
    qn = q(n,:);
    A = q2A(qn);
    % Body axes
    B = rot * A * rot';
    % Position
    pos = y(n,1:3);
    radial = pos / norm(pos);
    tangential = y(n,4:6) / norm(y(n,4:6));
    normal = cross(radial,tangential);
    RTN = [radial' tangential' normal'];
end

%% Problem 2
% Initial conditions
perturbation = 0.01;
w0x = [1; perturbation; perturbation];
w0y = [perturbation; 1; perturbation];
w0z = [perturbation; perturbation; 1];

w0Mat = {w0x, w0y, w0z};
q0 = [0; 0; 0; 1];

for n = 1:3
    w0 = w0Mat{n};
    
    [q,w] = kinQuaternionRK4(q0,w0,Ix,Iy,Iz,tFinal,tStep);

    figure(3)
    subplot(2,1,1)
    plot(t, w)
    xlabel('Time [s]')
    ylabel('Angular velocity (rad/s)')
    legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
        'Location','Southeast')

    subplot(2,1,2)
    plot(t, q)
    xlabel('Time [s]')
    ylabel('Quaternion')
    legend('q_{1}','q_{2}','q_{3}','q_{4}', ...
        'Location','Southeast')
    saveas(3, ['Images/ps4_problem2a_', sprintf('%i',n), '.png'])
end

%% Problem 3(c)
eulerAngle0 = [0; 360; 0];
w0 = [0.1; 0.1; 0.1; 100];
state0 = [eulerAngle0;w0];

tFinal = 60;
tStep = 0.1;
t = 0:tStep:tFinal;

M = [0; 0; 0; 0];
r = [0; 0; 1];
Ir = 100;

tspan = 0:tStep:tFinal;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinWheel(t,state,M,r,Ix,Iy,Iz,Ir), ...
    tspan,state0,options);

eulerAngle = wrapTo360(rad2deg(state(:,1:3)));

figure()
plot(t,state(:,4:6),'LineWidth',1)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(1,'Images/ps4_problem3c_velocity.png')

figure()
plot(t,state(:,1:3),'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel(['Euler Angle [' char(176) ']'])
saveas(1,'Images/ps4_problem3c_angle.png')

%% Problem 3(d)
eulerAngle0 = [0; 360; 0];
w0 = [0.1; 0.1; 0.1; 100];
state0 = [eulerAngle0;w0];

tFinal = 120;
tStep = 0.1;
t = 0:tStep:tFinal;

M = [0; 0; 0; 0];
r = [0; 1; 0];
Ir = 100;

tspan = 0:tStep:tFinal;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinWheel(t,state,M,r,Ix,Iy,Iz,Ir), ...
    tspan,state0,options);

eulerAngle = wrapTo360(rad2deg(state(:,1:3)));

figure()
plot(t,state(:,4:6),'LineWidth',1)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(1,'Images/ps4_problem3d_velocity.png')

figure()
plot(t,state(:,1:3),'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel(['Euler Angle [' char(176) ']'])
saveas(1,'Images/ps4_problem3d_angle.png')

%% Problem 3(e)
eulerAngle0 = [0; 360; 0];
w0 = [0.1; 0.1; 0.1; 100];
state0 = [eulerAngle0;w0];

tFinal = 120;
tStep = 0.1;
t = 0:tStep:tFinal;

M = [0; 0; 0; 0];
r = [sqrt(2)/2; sqrt(2)/2; 0];
Ir = 100;

tspan = 0:tStep:tFinal;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinWheel(t,state,M,r,Ix,Iy,Iz,Ir), ...
    tspan,state0,options);

eulerAngle = wrapTo360(rad2deg(state(:,1:3)));

figure()
plot(t,state(:,4:6),'LineWidth',1)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(1,'Images/ps4_problem3e_velocity.png')

figure()
plot(t,state(:,1:3),'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel(['Euler Angle [' char(176) ']'])
saveas(1,'Images/ps4_problem3e_angle.png')

%% Problem 4
% Get key orbital parameters
mu_E = 3.986e5; %km^3/s^2
P = 2*pi*sqrt(a^3/mu_E); %s
n_mm = 2*pi/P; %rad/s

r = y(:,1:3);
