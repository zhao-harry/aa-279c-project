%% Magnetic Field Investigation
close all; clear; clc

%% Playground
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

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
m = m_max*m_direction/norm(m_direction); % Arbitrary sat dipole
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

% figure()
% plot(t / 3600,state(:,7:9))
% xlabel('Time [h]')
% ylabel('Angular Velocity in Principal Axes [rad/s]')
% legend('\omega_{x}','\omega_{y}','\omega_{z}')
% 
% figure()
% plot(t / 3600,Mgg)
% xlabel('Time [h]')
% ylabel('Gravity Gradient Torque in Principal Axes [Nm]')
% legend('M_{x}','M_{y}','M_{z}')
% 
% figure()
% plot(t / 3600,Md)
% xlabel('Time [h]')
% ylabel('Drag Torque in Principal Axes [Nm]')
% legend('M_{x}','M_{y}','M_{z}')
% 
% figure()
% plot(t / 3600,Msrp)
% xlabel('Time [h]')
% ylabel('Solar Radiation Pressure Torque in Principal Axes [Nm]')
% legend('M_{x}','M_{y}','M_{z}')
% 
% figure()
% plot(t / 3600,Mm)
% xlabel('Time [h]')
% ylabel('Magnetic Field Torque in Principal Axes [Nm]')
% legend('M_{x}','M_{y}','M_{z}')

B_vec = zeros(size(state(:,1:3)));
for i = 1:length(t)
    r = state(i,1:3)';
    v = state(i,1:3)';
    radial = r / norm(r);
    rEarth = state(i,13:15)';
    A_ECI2P = e2A(state(i,10:12));
    
    M = magFieldTorque(m,r,state(10:12),t(i),6378.1,UT1);
    Mm(i,1:3) = M;
end

R = r0;
rE = rEarth;
rE3_B0 = 7.943e15; % Wb * m
GMST = time2GMST(0,UT12MJD(UT1));
[lat,lon,~] = ECEF2Geoc(ECI2ECEF(R,GMST),t);
lambda = lat;
phi = lon;
theta_m = pi/2 - phi;

dalphadt = 7.292115827689e-5;
t0 = 0;

alpha_m = GMST + dalphadt*t0 + lambda;

mE = [sin(theta_m)*cos(alpha_m);
         sin(theta_m)*sin(alpha_m);
         cos(theta_m)];

B1 = magFieldEarthDipole(R, rE, rE3_B0);
B1_norm = norm(B1);
[B2_1, B2_2, B2_3] = magFieldEarth(R, lambda, theta_m, norm(rE));
B2 = [B2_1; B2_2; B2_3];
B2_norm = norm(B2);

Mm_test1 = cross(m, B2);

tFinal = 60000;
tStep = 10;
t = 0:tStep:tFinal;

B_ECEF = zeros(size(state(:,1:3)));
% F = zeros(size(state(:,1:3)));
F = zeros(size(state(:,1)));

for n = 1:length(t)
    r = state(i,1:3)';
    v = state(i,1:3)';
    radial = r / norm(r);
    rEarth = state(i,13:15)';
    A_ECI2P = e2A(state(i,10:12));

    [M, B, F_val, XYZ] = magFieldTorque(m,r,state(n,10:12),t(n),6378.1,UT1);
    B_ECEF(n,1:3) = B;
    F(n, :) = F_val*1e-9;
end

t_days = t/86400;
figure(1)
hold on
% plot(t_days, B_ECEF)
% plot(t_days, F)
% legend(["B_x", "B_y", "B_z", "F_x", "F_y", "F_z"])
plot(t_days, vecnorm(B_ECEF'))
plot(t_days, F)
% plot(t_days, vecnorm(F'))
% legend("B")
hold off
