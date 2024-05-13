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
IR = Ix;
IT = Iz;
IN = Iy;

kT = (IN - IR) / IT;
kR = (IN - IT) / IR;

plotGravGradStability(kR,kT,'Nominal','Images/ps5_problem1a.png');

%% Problem 1(b) (Unstable, Unperturbed)
tFinal = 6000 * 5; % 5 orbits
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
A_Nominal = [-radial -normal -tangential]';

state0 = zeros(12,1);
state0(1:6) = y;
state0(7:9) = [0; -n; 0];
state0(10:12) = A2e(A_Nominal);

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
    saveas(gcf,'Images/ps5_problem1b_angvel_unperturbed.png')
end

figure()
plot(t / 3600,wrapToPi(state(:,10:12)))
xlabel('Time [h]')
ylabel('Euler Angles in Principal Axes [rad]')
legend('\phi','\theta','\psi')
if savePlot == true
    saveas(gcf,'Images/ps5_problem1b_angle_unperturbed.png')
end

%% Problem 1(b) (Unstable, Perturbed)
tFinal = 6000 * 3; % 2 orbits
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
A_Nominal = [-radial -normal -tangential]';

state0 = zeros(12,1);
state0(1:6) = y;
state0(7:9) = [0; -n; 0] * 1.01;
state0(10:12) = A2e(A_Nominal) + pi * [0.01; 0.01; 0.01];

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
plot(t / 3600,wrapToPi(state(:,10:12)))
xlabel('Time [h]')
ylabel('Euler Angles in Principal Axes [rad]')
legend('\phi','\theta','\psi')
if savePlot == true
    saveas(gcf,'Images/ps5_problem1b_angle.png')
end

%% Problem 1(c)
IR = Ix;
IT = Iy;
IN = Iz;
kT = (IN - IR) / IT;
kR = (IN - IT) / IR;

plotGravGradStability(kR,kT, ...
    'Principal XYZ aligned with RTN', ...
    'Images/ps5_problem1c.png');

%% Problem 1(c) (Stable, Perturbed)
tFinal = 6000 * 10; % 10 orbits
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
    saveas(gcf,'Images/ps5_problem1c_angvel.png')
end

figure()
plot(t / 3600,wrapToPi(state(:,10:12)))
xlabel('Time [h]')
ylabel('Euler Angles in Principal Axes [rad]')
legend('\phi','\theta','\psi')
if savePlot == true
    saveas(gcf,'Images/ps5_problem1c_angle.png')
end

%% Problem 3
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
muE = 3.986 * 10^5; % km^3 / s^2
n = sqrt(muE / a^3);

% Compute initial position and attitude
y = oe2eci(a,e,i,O,w,nu);
r0 = y(1:3);
v0 = y(4:6);
h = cross(r0,v0);
radial = r0 / norm(r0);
normal = h / norm(h);
tangential = cross(normal,radial);
A_RTN = [radial tangential normal]';

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
state0(7:9) = [0; 0; n];
state0(10:12) = A2e(A_RTN);
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
P = 1358/3E8;
S_sat = 24.92;
m_max = 4*pi*1e-7 * S_sat * 0.1;
m_direction_body = [1; 0; 0];
m_direction = rot * m_direction_body;
m = m_max*m_direction/norm(m_direction); % Arbitrarily defined satellite dipole for now
UT1 = [2024 1 1];

% Run numerical method
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) orbitTorque(t,state,Ix,Iy,Iz, ...
    CD,Cd,Cs,P,m,UT1, ...
    barycenter,normal,area,cmP,n), ...
    tspan,state0,options);

% Compute torques (since ode113 does not allow returning these)
Rx = [1 0 0; 0 cosd(23.5) -sind(23.5); 0 sind(23.5) cosd(23.5)];
c = zeros(size(state(:,1:3)));
Mgg = zeros(size(state(:,1:3)));
Md = zeros(size(state(:,1:3)));
Msrp = zeros(size(state(:,1:3)));
Mm = zeros(size(state(:,1:3)));

for i = 1:length(t)
    r = state(i,1:3)';
    v = state(i,1:3)';
    radial = r / norm(r);
    rEarth = state(i,13:15)';
    A_ECI2P = e2A(state(i,10:12));

    c(i,1:3) = A_ECI2P * radial;
    Mgg(i,1:3) = gravGradTorque(Ix,Iy,Iz,n,c(i,1:3));
    
    [~,density] = atmosnrlmsise00(1000 * (norm(r) - 6378.1),0,0,2000,1,0);
    rho = density(6);
    vPrincipal = A_ECI2P * (v + cross([0; 0; 7.2921159E-5],r));
    [~,M] = drag(vPrincipal,rho,CD,barycenter,normal,area,cmP);
    Md(i,1:3) = M;

    s = A_ECI2P * (-Rx * rEarth - r);
    [~,M] = srp(s,P,Cd,Cs,barycenter,normal,area,cmP);
    Msrp(i,1:3) = M;

    M = magFieldTorque(m,r,state(i,10:12),t(i),6378.1,UT1);
    Mm(i,1:3) = M;
end

figure()
plot(t / 3600,state(:,7:9))
xlabel('Time [h]')
ylabel('Angular Velocity in Principal Axes [rad/s]')
legend('\omega_{x}','\omega_{y}','\omega_{z}')
if savePlot == true
    saveas(gcf,'Images/ps5_problem3_angvel.png')
end

figure()
plot(t / 3600,Mgg)
xlabel('Time [h]')
ylabel('Gravity Gradient Torque in Principal Axes [Nm]')
legend('M_{x}','M_{y}','M_{z}')
if savePlot == true
    saveas(gcf,'Images/ps5_problem3_grav.png')
end

figure()
plot(t / 3600,Md)
xlabel('Time [h]')
ylabel('Drag Torque in Principal Axes [Nm]')
legend('M_{x}','M_{y}','M_{z}')
if savePlot == true
    saveas(gcf,'Images/ps5_problem3_drag.png')
end

figure()
plot(t / 3600,Msrp)
xlabel('Time [h]')
ylabel('Solar Radiation Pressure Torque in Principal Axes [Nm]')
legend('M_{x}','M_{y}','M_{z}')
if savePlot == true
    saveas(gcf,'Images/ps5_problem3_srp.png')
end

figure()
plot(t / 3600,Mm)
xlabel('Time [h]')
ylabel('Magnetic Field Torque in Principal Axes [Nm]')
legend('M_{x}','M_{y}','M_{z}')
if savePlot == true
    saveas(gcf,'Images/ps5_problem3_mag.png')
end

%%
figure()
plot(t, state(:,10:12))

%% Problem 3 Maximum Torques
% Parameters
CD = 2;
Cd = 0; Cs = 0.9;
q = Cd + Cs;
P = 1358/3E8;
S_sat = 24.92;
m_max = 4*pi*1e-7 * S_sat * 0.1;
UT1 = [2024 1 1];
rE3_B0 = 7.943e15; %Wb*km

r_norm = norm(state(1,1:3));
Mgg_max = 3/2 * muE/(r_norm^3) * abs(max([Ix Iy Iz]) - min([Ix Iy Iz]));
Mm_max = 2*m_max*rE3_B0/((r_norm*1e3)^3);
Msrp_max = 0;
Md_max = 0;

vMax = max(vecnorm(state(:,4:6)')) * 1e3;

for n = 1:length(area)
    Msrp_max = Msrp_max + P*area(n)*(1+q) * norm(barycenter(:,n) - cmP);
    Md_max = Md_max + 0.5*rho*CD*vMax^2*area(n) * norm(barycenter(:,n) - cmP);
end

fprintf("Maximum expected values: \n" + ...
        "M_gg: %d Nm \n" + ...
        "M_m: %d Nm \n" + ...
        "M_srp: %d Nm \n" + ...
        "M_d: %d Nm\n", Mgg_max, Mm_max, Msrp_max, Md_max);