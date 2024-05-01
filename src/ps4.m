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
    saveas(gcf,'Images/ps4_problem1a_angvel.png')
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
momentumPlot = 'Images/ps4_problem3c_angmom.png';
velocityPlot = 'Images/ps4_problem3c_angvel.png';
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
momentumPlot = 'Images/ps4_problem3d_angmom.png';
velocityPlot = 'Images/ps4_problem3d_angvel.png';
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
momentumPlot = 'Images/ps4_problem3e_angmom.png';
velocityPlot = 'Images/ps4_problem3e_angvel.png';
anglePlot = 'Images/ps4_problem3e_angle.png';
plotPS4Problem3(eulerAngle0,w0,tStep,tFinal, ...
                M,r,Ix,Iy,Iz,Ir, ...    
                momentumPlot,velocityPlot,anglePlot,savePlot);

%% Problem 4
tFinal = 60;
tStep = 0.1;
tspan = 0:tStep:tFinal;

a = 7125.48662; % km
e = 0;
i = 98.40508; % degree
O = -19.61601; % degree
w = 89.99764; % degree
nu = -89.99818; % degree
muE = 3.986e5;

n = sqrt(muE / a^3);

y = oe2eci(a,e,i,O,w,nu);
r0 = y(1:3);
v0 = y(4:6);
h = cross(r0,v0);
radial = r0 / norm(r0);
normal = h / norm(h);
tangential = cross(normal,radial);
A_RTN = [radial tangential normal];

state0 = zeros(12,1);
state0(1:6) = y;
state0(7:9) = [0; 0; n];
state0(10:12) = A2e(A_RTN);

options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) gravGrad(t,state,Ix,Iy,Iz,n), ...
        tspan,state0,options);

c = zeros(size(state(:,1:3)));
M = zeros(size(state(:,1:3)));
for i = 1:length(t)
    r = state(i,1:3);
    v = state(i,4:6);
    h = cross(r,v);
    radial = r / norm(r);
    normal = h / norm(h);
    tangential = cross(normal,radial);
    A_RTN = [radial' tangential' normal'];
    A_ECI2P = e2A(state(i,10:12));
    A = A_ECI2P * A_RTN';
    c(i,1:3) = A(:,1);
    M(i,1:3) = gravGradTorque(Ix,Iy,Iz,n,c(i,1:3));
end

% plot(t,c)

% plot(t,M)

% plot(t,state(:,7:9))
% xlabel('Time [s]')
% ylabel('Angular Velocity in Principal Axes [rad/s]')
% legend('\omega_{x}','\omega_{y}','\omega_{z}')

%%
stateT = state(end,:);
r = stateT(1:3);
v = stateT(4:6);
h = cross(r,v);
radial = r / norm(r);
normal = h / norm(h);
tangential = cross(normal,radial);
A_RTN = [radial' tangential' normal']
A_ECI2P = e2A(stateT(10:12))

%% Problem 4(c)
% M = 3 * muE / a^3 * [(Iz - Iy) * c(2) * c(3); ...
%                     (Ix - Iz) * c(3) * c(1); ...
%                     (Iy - Ix) * c(1) * c(2)]
