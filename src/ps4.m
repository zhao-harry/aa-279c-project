close all; clear; clc
savePlot = false;

%% Import mass properties
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

%% Problem 1(a)
% Find Euler kinematics
tFinal = 60;
tStep = 0.01;
tspan = 0:tStep:tFinal;

eulerAngle0 = [0; 0; 0];
w0 = [0; 0; 1];
state0 = [eulerAngle0;w0];

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
if savePlot
    saveas(gcf,'Images/ps4_problem1a_velocity.png')
end

figure()
plot(t,state(:,1:3),'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel(['Euler Angle [' char(176) ']'])
if savePlot
    saveas(gcf,'Images/ps4_problem1a_angle.png')
end

%% Problem 1(b)
a = 7125.48662; % km
e = 0.0011650;
i = 98.40508; % degree
O = -19.61601; % degree
w_deg = 89.99764; % degree
nu = -89.99818; % degree

tFinal = 6000;
tStep = 0.1;
tspan = 0:tStep:tFinal;
tTrunc = 300;
nTrunc = find((tspan == tTrunc) == 1);

[~,y] = plotECI(a,e,i,O,w_deg,nu,tspan);
close all
format long

% Initialize w0 to be aligned with normal
r0 = y(1,1:3);
v0 = y(1,4:6);
h = cross(r0,v0);
radial = r0 / norm(r0);
normal = h / norm(h);
tangential = cross(normal,radial);
A_RTN = [radial' tangential' normal'];
w0_RTN = [0; 0; 0.1];
euler0_RTN = A2e(A_RTN);
state0 = [euler0_RTN; w0_RTN];
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinEulerAngle(t,state,Ix,Iy,Iz), ...
    tspan,state0,options);

w_RTN = nan(size(state(:,1:3)));
euler_RTN = nan(size(state(:,1:3)));

for n = 1:length(t)
    pos = y(n,1:3);
    vel = y(n,4:6);
    h = cross(pos,vel);
    radial = pos / norm(pos);
    normal = h / norm(h);
    tangential = cross(normal,radial);
    A_RTN = [radial' tangential' normal'];

    % Get rotation matrixes (to ECI)
    euler = state(n,1:3);
    w_principal = state(n,4:6)';
    A_principal = e2A(euler);

    A_P2R = A_RTN * A_principal';

    w_RTN(n,:) = A_P2R*w_principal;
    euler_RTN(n,:) = A2e(A_P2R');
end


figure(2)
plot(t, w_RTN)
legend("\omega_{R}", "\omega_{T}", "\omega_{N}")
if savePlot == true
    saveas(2, 'Images/ps4_problem1b_angvel.png')
end


figure(3)
plot(t(1:nTrunc), euler_RTN(1:nTrunc,:))
title("Euler Angles from principal to RTN frame")
legend("\phi", "\theta", "\psi")
if savePlot == true
    saveas(3, 'Images/ps4_problem1b_euler.png')
end

%% Problem 2
% Initial conditions
perturbation = 0.001;
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

    figure()
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
    if savePlot
        saveas(gcf, ['Images/ps4_problem2a_', sprintf('%i',n), '.png'])
    end
end

%% Problem 3(c)
eulerAngle0 = [0; 0; 0];
w0 = [0.01; 0.01; 0.01; 100];
M = [0; 0; 0; 0];
r = [0; 0; 1];
Ir = 100;
tFinal = 60;
tStep = 0.1;
momentumPlot = 'Images/ps4_problem3c_momentum.png';
velocityPlot = 'Images/ps4_problem3c_velocity.png';
anglePlot = 'Images/ps4_problem3c_angle.png';
plotPS4Problem3(eulerAngle0,w0,tStep,tFinal, ...
                M,r,Ix,Iy,Iz,Ir, ...    
                momentumPlot,velocityPlot,anglePlot,savePlot);

%% Problem 3(d)
eulerAngle0 = [0; 0; 0];
w0 = [0.01; 0.25; 0.01; 100];
M = [0; 0; 0; 0];
r = [0; 1; 0];
Ir = 100;
tFinal = 120;
tStep = 0.1;
momentumPlot = 'Images/ps4_problem3d_momentum.png';
velocityPlot = 'Images/ps4_problem3d_velocity.png';
anglePlot = 'Images/ps4_problem3d_angle.png';
plotPS4Problem3(eulerAngle0,w0,tStep,tFinal, ...
                M,r,Ix,Iy,Iz,Ir, ...    
                momentumPlot,velocityPlot,anglePlot,savePlot);

%% Problem 3(e)
eulerAngle0 = [0; pi/2; 0];
w0 = [0.1; 0.1; 0.1; 100];
M = [0; 0; 0; 0];
r = [sqrt(2)/2; sqrt(2)/2; 0];
Ir = 100;
tFinal = 120;
tStep = 0.1;
momentumPlot = 'Images/ps4_problem3e_momentum.png';
velocityPlot = 'Images/ps4_problem3e_velocity.png';
anglePlot = 'Images/ps4_problem3e_angle.png';
plotPS4Problem3(eulerAngle0,w0,tStep,tFinal, ...
                M,r,Ix,Iy,Iz,Ir, ...    
                momentumPlot,velocityPlot,anglePlot,savePlot);

%% Problem 4
% Get key orbital parameters
orbitParams.a = 7125.48662; % km
orbitParams.e = 0.0011650;
orbitParams.i = 98.40508; % degree
orbitParams.O = -19.61601; % degree
orbitParams.w = 89.99764; % degree
orbitParams.nu = -89.99818; % degree

% Get data
tFinal = 600;
tStep = 1;
t = 0:tStep:tFinal;
w0 = w0_RTN; %from 1b
eulerAng0 = euler0_RTN; % from 1b

[w, eulerAngs, M_gg, y] = eulerEquationGravGradRK4(w0,eulerAng0,Ix,Iy,Iz,orbitParams,tFinal,tStep);
eulerAngs = wrapTo360(rad2deg(eulerAngs));

figure(5)
plot(t, w)