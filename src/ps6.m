close all; clear; clc

% From PS6 onward, we use Simulink to model the spacecraft
savePlots = false;

%% Import mass properties
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

%% Problem 2-3
tFinal = 6000;
tStep = 1;
tspan = 0:tStep:tFinal;

% Satellite orbit initial conditions
a = 7125.48662; % km
e = 0;
i = 98.40508; % degree
O = -19.61601; % degree
w = 89.99764; % degree
nu = -89.99818; % degree
muE = 3.986e5; % km^3 / s^2
n = sqrt(muE / a^3);

% Compute initial position and attitude
y = oe2eci(a,e,i,O,w,nu);
r0 = y(1:3);
v0 = y(4:6);
h = cross(r0,v0);
radial = r0 / norm(r0);
normal = h / norm(h);
tangential = cross(normal,radial);
A_Nominal = [-radial -normal -tangential]';

% Earth orbit initial conditions
aE = 149.60E6; % km
eE = 0.0167086;
iE = 7.155; % degree
OE = 174.9; % degree
wE = 288.1; % degree
nuE = 0;
muSun = 1.327E11; % km^3 / s^2
nE = sqrt(muSun / aE^3);
ySun = oe2eci(aE,eE,iE,OE,wE,nuE);

% Initial conditions
state0 = zeros(12,1);
state0(1:6) = y;
state0(7:9) = [0; -n; 0];
state0(10:12) = A2e(A_Nominal);
state0(13:18) = ySun;

% Properties
[barycenter,normal,area] = surfaces('res/area.csv',rot');
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);
[rot,~] = eig(I);
cmP = rot' * cm;

% Parameters
CD = 2;
Cd = 0; Cs = 0.9;
P = 1358 / 3e8;
S_sat = 24.92;
m_max = 4e-7 * pi * S_sat * 0.1;
m_direction_body = [1; 0; 0];
m_direction = rot * m_direction_body;
m = m_max * m_direction / norm(m_direction); % Arbitrary sat dipole
UT1 = [2024 1 1];

% Run numerical method
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) orbitTorque(t,state,Ix,Iy,Iz, ...
    CD,Cd,Cs,P,m,UT1, ...
    barycenter,normal,area,cmP,n), ...
    tspan,state0,options);

eulerError = zeros(size(state(:,10:12)));
eulerTarget = zeros(size(state(:,10:12)));
for i = 1:length(t)
    r = state(i,1:3)';
    v = state(i,4:6)';
    h = cross(r,v);
    radial = r / norm(r);
    normal = h / norm(h);
    tangential = cross(normal,radial);
    A_Target = [-radial -normal -tangential]'; % ECI -> RTN
    A_Sat = e2A(state(i,10:12)); % ECI -> Principal
    eulerTarget(i,1:3) = A2e(A_Target);
    eulerError(i,1:3) = A2e(A_Sat * A_Target');
end

figure()
plot(t / 3600,wrapToPi(state(:,10:12)))
xlabel('Time [h]')
ylabel('Euler Angle (Principal) [rad]')
legend('\phi','\theta','\psi')
saveAsBool(gcf,'Images/ps6_problem2_principal.png',savePlots)

figure()
plot(t / 3600,wrapToPi(eulerTarget))
xlabel('Time [h]')
ylabel('Euler Angles (Target) [rad]')
legend('\phi','\theta','\psi')
saveAsBool(gcf,'Images/ps6_problem2_target.png',savePlots)

figure()
plot(t / 3600,wrapToPi(eulerError))
xlabel('Time [h]')
ylabel('Euler Angle Error [rad]')
legend('\phi','\theta','\psi')
saveAsBool(gcf,'Images/ps6_problem2_error.png',savePlots)

%% Problem 5



numReadings = 2;

v = rand([3, numReadings]);
v  = v ./ vecnorm(v);

m = nan(size(v));

for n = 1:numReadings
    m(:,n) = A_Nominal * v(:,n);
end

qReal = A2q(A_Nominal);

w = [100 100 ones([1, numReadings-2])];

m1 = m(:,1);
m2 = m(:,2);
v1 = v(:,1);
v2 = v(:,2);

A1 = DAD_twoVecs(m1, m2, v1, v2);
q1 = A2q(A1)

A2 = DAD(m, v);
q2 = A2q(A2)

q3 = qMethod(m,v,w);