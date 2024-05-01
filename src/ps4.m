close all; clear; clc
savePlot = true;

%% Import mass properties
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

%% Problem 1(a)
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

% Plot
figure()
plot(t,state(:,4:6),'LineWidth',1)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular Velocity (\omega) [' char(176) '/s]'])
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
e = 0;
i = 98.40508; % degree
O = -19.61601; % degree
w_deg = 89.99764; % degree
nu = -89.99818; % degree

muE = 3.986 * 10^5;

n = sqrt(muE / a^3);

tFinal = 6000;
tStep = 0.1;
tspan = 0:tStep:tFinal;
tTrunc = 300;
nTrunc = find((tspan == tTrunc) == 1);

[~,y] = plotECI(a,e,i,O,w_deg,nu,tspan);
close all
format long

% Initialize angular velocity aligned with normal
r0 = y(1,1:3);
v0 = y(1,4:6);
h = cross(r0,v0);
radial = r0 / norm(r0);
normal = h / norm(h);
tangential = cross(normal,radial);
A_RTN = [radial' tangential' normal']';
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
    tangential = tangential / norm(tangential);
    A_RTN = [radial' tangential' normal']';

    % Get rotation matrixes (to ECI)
    euler = state(n,1:3);
    w_principal = state(n,4:6)';
    A_principal = e2A(euler);
    A_P2R = A_RTN * A_principal';
    w_RTN(n,:) = A_P2R*w_principal;
    euler_RTN(n,:) = A2e(A_P2R');
end

figure()
plot(t, rad2deg(w_RTN))
xlabel('Time [s]')
ylabel(['Angular Velocity, RTN Frame [' char(176) '/s]'])
legend('\omega_{R}','\omega_{T}','\omega_{N}')
if savePlot == true
    saveas(gcf,'Images/ps4_problem1b_angvel.png')
end

figure()
plot(t(1:nTrunc), wrapTo180(rad2deg(euler_RTN(1:nTrunc,:))))
xlabel('Time [s]')
ylabel(['Euler Angle, RTN Frame [' char(176) ']'])
legend('\phi','\theta','\psi')
if savePlot == true
    saveas(gcf,'Images/ps4_problem1b_euler.png')
end

%% Problem 2
% Initial conditions
perturbation = 0.001;
w0x = [1; perturbation; perturbation];
w0y = [perturbation; 1; perturbation];
w0z = [perturbation; perturbation; 1];

w0Mat = {w0x, w0y, w0z};
eulerAngle0 = [0; 0; 0];
tStep = 0.01;
tFinal = 60;

for n = 1:3
    w0 = w0Mat{n};
    state0 = [eulerAngle0;w0];

    tspan = 0:tStep:tFinal;
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t,state] = ode113(@(t,state) kinEulerAngle(t,state,Ix,Iy,Iz), ...
        tspan,state0,options);

    eulerAngle = wrapTo180(rad2deg(state(:,1:3)));
    w = wrapTo180(rad2deg(state(:,4:6)));

    figure()
    subplot(2,1,1)
    plot(t,w)
    xlabel('Time [s]')
    ylabel(['Angular Velocity [' char(176) '/s]'])
    legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
        'Location','Southeast')

    subplot(2,1,2)
    plot(t,eulerAngle)
    xlabel('Time [s]')
    ylabel(['Euler Angle [' char(176) ']'])
    legend('\phi','\theta','\psi', 'Location','Southeast')
    if savePlot
        saveas(gcf,['Images/ps4_problem2a_' sprintf('%i',n) '.png'])
    end
end

%% Momentum wheel setup
% Based on RSI 68
mr = 8.9; % kg
r = 0.347 / 2; % m
Ir = mr * r^2;
wrRPM = 2500; % RPM
wr = wrRPM * 0.1047198;

% No external torques
M = [0; 0; 0; 0];

% Time
tFinal = 600;
tStep = 0.1;

%% Problem 3(b)
eulerAngle0 = [0; 0; 0];
w0 = [0; 0; 0; wr];
r = [0; 0; 1];

momentumPlot = 'Images/ps4_problem3b_angmom.png';
velocityPlot = 'Images/ps4_problem3b_angvel.png';
anglePlot = 'Images/ps4_problem3b_angle.png';
plotPS4Problem3(eulerAngle0,w0,tStep,tFinal, ...
                M,r,Ix,Iy,Iz,Ir, ...    
                momentumPlot,velocityPlot,anglePlot,savePlot);

%% Problem 3(c)
eulerAngle0 = [0; 0; 0];
w0 = [0.01; 0.01; 0.01; wr];
r = [0; 0; 1];

momentumPlot = 'Images/ps4_problem3c_angmom.png';
velocityPlot = 'Images/ps4_problem3c_angvel.png';
anglePlot = 'Images/ps4_problem3c_angle.png';
plotPS4Problem3(eulerAngle0,w0,tStep,tFinal, ...
                M,r,Ix,Iy,Iz,Ir, ...    
                momentumPlot,velocityPlot,anglePlot,savePlot);

%% Problem 3(d)
eulerAngle0 = [0; 0; 0];
w0 = [0.01; 0.25; 0.01; wr];
r = [0; 1; 0];

momentumPlot = 'Images/ps4_problem3d_angmom.png';
velocityPlot = 'Images/ps4_problem3d_angvel.png';
anglePlot = 'Images/ps4_problem3d_angle.png';
plotPS4Problem3(eulerAngle0,w0,tStep,tFinal, ...
                M,r,Ix,Iy,Iz,Ir, ...    
                momentumPlot,velocityPlot,anglePlot,savePlot);

%% Problem 3(e)
eulerAngle0 = [0; 0; 0];
w0 = [0.1; 0.1; 0.1; wr];
r = [sqrt(2)/2; sqrt(2)/2; 0];

momentumPlot = 'Images/ps4_problem3e_angmom.png';
velocityPlot = 'Images/ps4_problem3e_angvel.png';
anglePlot = 'Images/ps4_problem3e_angle.png';
plotPS4Problem3(eulerAngle0,w0,tStep,tFinal, ...
                M,r,Ix,Iy,Iz,Ir, ...    
                momentumPlot,velocityPlot,anglePlot,savePlot);

%% Problem 4
tFinal = 6000;
tStep = 1;
tspan = 0:tStep:tFinal;

a = 7125.48662; % km
e = 0;
i = 98.40508; % degree
O = -19.61601; % degree
w = 89.99764; % degree
nu = -89.99818; % degree
muE = 3.986 * 10^5;
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
state0(10:12) = A2e(A_RTN');

options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) gravGrad(t,state,Ix,Iy,Iz,n), ...
        tspan,state0,options);

c = zeros(size(state(:,1:3)));
M = zeros(size(state(:,1:3)));
for i = 1:length(t)
    r = state(i,1:3);
    radial = r / norm(r);
    A_ECI2P = e2A(state(i,10:12));
    c(i,1:3) = A_ECI2P * radial';
    M(i,1:3) = gravGradTorque(Ix,Iy,Iz,n,c(i,1:3));
end

% plot(t,M)
% xlabel('Time [s]')
% ylabel('Torque in Principal Axes [rad/s]')
% legend('M_{x}','M_{y}','M_{z}')
% ylim([-1e-5 1e-5])

plot(t,state(:,7:9))
xlabel('Time [s]')
ylabel('Angular Velocity in Principal Axes [rad/s]')
legend('\omega_{x}','\omega_{y}','\omega_{z}')

%%
figure()
a = 7125.48662; % km
e = 0;
i = 98.40508; % degree
O = -19.61601; % degree
w = 89.99764; % degree
nu = -89.99818; % degree
muE = 3.986 * 10^5;
y = state(:,1:6);
plotECI(a,e,i,O,w,nu,tspan);
hold on
tLen = length(t);
for i = 1:50:tLen
    r = y(i,1:3);
    v = y(i,4:6);
    radial = r / norm(r);
    h = cross(r,v);
    normal = h / norm(h);
    tangential = cross(normal,radial);
    RTN = [radial' tangential' normal'];
    A_P = e2A(state(i,10:12))';
    plotTriad(gca,r,A_P,500,'b');
end
hold off
