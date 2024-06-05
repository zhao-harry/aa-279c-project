%% Simulink Variables 
% Values needed to run before simulink model works!

%% Settings
% What kind of measurement method?
% measType = "dad";
% measType = "q";
% measType = "kin";
measType = "MEKF";

% If dad, do we use fictitious measurements?
useFict = true;
% useFict = false;

% Do we use linear or nonlinear control model?
% useLinearModel = true;
useLinearModel = false;

%% General Constants
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

constants = struct();
constants.Ix = Ix; constants.Iy = Iy; constants.Iz = Iz;
constants.A_Body2P = rot;
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
P0 = eye(6);

%% Sensor information
num_stars = 10;
sensor_weights = [50 1]; %[starTracker, sunSensor]
sun_sensor_error = deg2rad(0.5);
star_tracker_error = deg2rad(0.01);
gyro_error = deg2rad(0.001);
gyro_bias = 5e-5;
star_tracker_normal_body = {[0; 1; 0], [0; -1; 0]};
star_tracker_FOV = deg2rad(20);
[~, indBest2Sensors] = maxk(sensor_weights, 2);
sensors_Q = eye(6)/100;
sensors_R = diag([repelem(star_tracker_error, 2*num_stars), sun_sensor_error]).^2;

sensors = struct();
sensors.weights = sensor_weights;
sensors.sunError = sun_sensor_error;
sensors.numStars = num_stars;
sensors.trackerError = star_tracker_error;
sensors.trackerFOV = star_tracker_FOV;
sensors.gyroError = gyro_error;
sensors.gyroBias = gyro_bias;
sensors.Q = sensors_Q;
sensors.R = sensors_R;
sensors_bus_info = Simulink.Bus.createObject(sensors);
sensors_bus = evalin('base', sensors_bus_info.busName);

%% Controller
% Actuators
IWheel = 0.119; % kg*m^2
LMaxWheel = 50; %N*m*s
dipole = [350; 350; 565]; %A*m^2
thrust = 1; % N
Isp = 220; %s

% Control parameters
f = 250;
Kp = [f^2/Ix;
         f^2/Iy;
         f^2/Iz];
Kd = [2*sqrt(Ix*(3*n^2*(Iz-Iy) + Kp(1)));
         2*sqrt(Iy*(3*n^2*(Ix-Iz) + Kp(2)));
         2*sqrt(Iz*(3*n^2*(Iy-Ix) + Kp(3)))];
A = 1/sqrt(3)*[-1 1 1 -1;
                      -1 -1 1 1;
                      1 1 1 1];
Lw0 = [0; 0; 0; 0];

% Actuator bus
control = struct();
control.IWheel = IWheel;
control.LMaxWheel = LMaxWheel;
control.magnetorquerDipole = dipole;
control.thrust = 1;
control.Isp = Isp;
control.Kp = Kp;
control.Kd = Kd;
control.A = A;
control.useLinearModel = useLinearModel;
control_bus_info = Simulink.Bus.createObject(control);
control_bus = evalin('base', control_bus_info.busName);

