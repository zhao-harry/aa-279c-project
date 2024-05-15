close all; clear; clc

%% Get variables for simulink
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

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

% Sensor information (assume 5 readings)
numReadings = 5;
sensor_weights = [100 100 ones([1, numReadings-2])];
[~, indBest2Sensors] = maxk(sensor_weights, 2);
ground_truth_vectors = rand([3, numReadings]);
ground_truth_vectors = ground_truth_vectors ./ norm(ground_truth_vectors);

% Get simulink vars
constants = struct();
constants.Ix = Ix; constants.Iy = Iy; constants.Iz = Iz;
constants.RE = 6378.1; %km
constants.n = n;
constants.cm = cmP;
constants.barycenter = barycenter;
constants.normal = normal;
constants.area = area;
constants.CD = CD;
constants.P = P;
constants.Cd = Cd;
constants.Cs = Cs;
constants.UT1 = UT1;
constants.m = m;
constants_bus_info = Simulink.Bus.createObject(constants);
constants_bus = evalin('base', constants_bus_info.busName);
rECI0 = r0; vECI0 = v0;
rSCI0 = ySun(1:3); vSCI0 = ySun(4:6);
q0 = A2q(A_Nominal);
w0 = [0, -n, 0];

sensors = struct();
sensors.weights = sensor_weights;
sensors.ground_truth_vectors = ground_truth_vectors;
sensors.indBest2Sensors = indBest2Sensors;
sensors_bus_info = Simulink.Bus.createObject(sensors);
sensors_bus = evalin('base', sensors_bus_info.busName);

% Settings
% measType = "dad";
measType = "q";
% measType = "kin";
useFict = true;
% useFict = false;


%% Plot

qVals = squeeze(out.q.data);
qMeasVals = squeeze(out.qMeasured.data);
timeVals = squeeze(out.q.time);

figure(1)
plot(timeVals, qVals)

figure(2)
plot(timeVals, qMeasVals-qVals)