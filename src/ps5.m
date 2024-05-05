close all; clear; clc
savePlot = true;

%% Import mass properties
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

%% Problem 1


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
[t,state] = ode113(@(t,state) gravGrad(t,state,Ix,Iy,Iz,n), ...
        tspan,state0,options);