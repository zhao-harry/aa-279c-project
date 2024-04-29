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

eulerAngle0 = [0; 1e-9; 0];
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
    % Get rotation matrixes (to ECI)
    euler = state(n,1:3);
    w_principal = state(n,4:6)';
    A_principal = euler2A(euler);

    pos = y(n,1:3);
    radial = pos / norm(pos);
    tangential = y(n,4:6) / norm(y(n,4:6));
    normal = cross(radial,tangential);
    A_RTN = [radial' tangential' normal'];

    A = A_RTN' * A_principal;
    w_RTN(n,:) = A*w_principal;
end

figure(2)
plot(t, w_RTN)

%% Problem 2
% Initial conditions
perturbation = 0.01;
w0x = [1; perturbation; perturbation];
w0y = [perturbation; 1; perturbation];
w0z = [perturbation; perturbation; 1];

w0Mat = {w0x, w0y, w0z};
eulerAngle0 = [0; 1e-9; 0];
tStep = 0.01;
tFinal = 30;

for n = 1:3
    w0 = w0Mat{n};
    state0 = [eulerAngle0;w0];
    
    tspan = 0:tStep:tFinal;
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t,state] = ode113(@(t,state) kinEulerAngle(t,state,Ix,Iy,Iz), ...
        tspan,state0,options);

    eulerAngle = wrapTo360(rad2deg(state(:,1:3)));
    w = state(:,4:6);

    figure(3)
    subplot(2,1,1)
    plot(t, w)
    xlabel('Time [s]')
    ylabel('Angular velocity (rad/s)')
    legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
        'Location','Southeast')

    subplot(2,1,2)
    plot(t, eulerAngle)
    xlabel('Time [s]')
    ylabel('Euler Angles')
    legend('\phi','\theta','\psi', 'Location','Southeast')
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
[t,state] = ode113(@(t,state) kinEulerAngleWheel(t,state,M,r,Ix,Iy,Iz,Ir), ...
    tspan,state0,options);

eulerAngle = wrapTo360(rad2deg(state(:,1:3)));

figure()
plot(t,state(:,4:6),'LineWidth',1)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(gcf,'Images/ps4_problem3c_velocity.png')

figure()
plot(t,state(:,1:3),'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel(['Euler Angle [' char(176) ']'])
saveas(gcf,'Images/ps4_problem3c_angle.png')

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
[t,state] = ode113(@(t,state) kinEulerAngleWheel(t,state,M,r,Ix,Iy,Iz,Ir), ...
    tspan,state0,options);

eulerAngle = wrapTo360(rad2deg(state(:,1:3)));

figure()
plot(t,state(:,4:6),'LineWidth',1)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(gcf,'Images/ps4_problem3d_velocity.png')

figure()
plot(t,state(:,1:3),'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel(['Euler Angle [' char(176) ']'])
saveas(gcf,'Images/ps4_problem3d_angle.png')

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
[t,state] = ode113(@(t,state) kinEulerAngleWheel(t,state,M,r,Ix,Iy,Iz,Ir), ...
    tspan,state0,options);

eulerAngle = wrapTo360(rad2deg(state(:,1:3)));

figure()
plot(t,state(:,4:6),'LineWidth',1)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(gcf,'Images/ps4_problem3e_velocity.png')

figure()
plot(t,state(:,1:3),'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel(['Euler Angle [' char(176) ']'])
saveas(gcf,'Images/ps4_problem3e_angle.png')

%% Problem 4
% Get key orbital parameters
mu_E = 3.986e5; %km^3/s^2
P = 2*pi*sqrt(a^3/mu_E); %s
n_mm = 2*pi/P; %rad/s

r = y(:,1:3);
