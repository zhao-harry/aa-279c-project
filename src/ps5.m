close all; clear; clc
savePlot = true;

%% Import mass properties
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);
IxB = I(1,1);
IyB = I(2,2);
IzB = I(3,3);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

%% Problem 1(a)
IR = Ix;
IT = Iy;
IN = Iz;
kT = (IN - IR) / IT;
kR = (IN - IT) / IR;

plotGravGradStability(kR,kT,'Principal','Images/ps5_problem1a.png');

%% Problem 1(b) (Stable, Perturbed)
tFinal = 6000 * 100; % 100 orbits
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
A_RTN = [radial tangential normal]';

state0 = zeros(12,1);
state0(1:6) = y;
state0(7:9) = [0; 0; n] * 1.01;
state0(10:12) = A2e(A_RTN) + pi * [0.01; 0.01; 0.01];

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

figure()
plot(t / 3600,state(:,7:9))
xlabel('Time [h]')
ylabel('Angular Velocity in Principal Axes [rad/s]')
legend('\omega_{x}','\omega_{y}','\omega_{z}')
if savePlot == true
    saveas(gcf,'Images/ps5_problem1b_angvel.png')
end

figure()
plot(t / 3600,state(:,10:12))
xlabel('Time [h]')
ylabel('Euler Angles in Principal Axes [rad]')
legend('\phi','\theta','\psi')
if savePlot == true
    saveas(gcf,'Images/ps5_problem1b_angle.png')
end

%% Problem 1(c)
IR = IzB;
IT = IyB;
IN = IxB;
kT = (IN - IR) / IT;
kR = (IN - IT) / IR;

plotGravGradStability(kR,kT,'Body','Images/ps5_problem1c.png');

%% Problem 3
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
A_RTN = [radial tangential normal]';

state0 = zeros(12,1);
state0(1:6) = y;
state0(7:9) = [0; 0; n];
state0(10:12) = A2e(A_Body);

options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) orbitTorque(t,state,Ix,Iy,Iz,n), ...
        tspan,state0,options);