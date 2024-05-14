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

%% Problem 2
% tFinal = 86400;
tFinal = 2000;
tStep = 5;
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
A_ECI2RTN = [radial tangential normal]';
A_RTN2ideal = [-1 0 0;
                        0 0 1;
                        0 1 0];
% A_RTN2ideal = -eye(3);
A_ideal0 = A_RTN2ideal*A_ECI2RTN;

% figure(1)
% hold on
% quiver3(0, 0, 0, 1, 0, 0, 'g')
% quiver3(0, 0, 0, 0, 1, 0, 'c')
% quiver3(0, 0, 0, 0, 0, 1, 'm')
% for i = 2%1:3
%     quiver3(0, 0, 0, A_ECI2RTN(i,1), A_ECI2RTN(i,2), A_ECI2RTN(i,3), 'r--')
%     quiver3(0, 0, 0, A_ideal0(i,1), A_ideal0(i,2), A_ideal0(i,3), 'b')
% end
% hold off

%%
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
stateOrb0(1:6) = y;
stateEul0(1:3) = [0; 0; n];
stateEul0(4:6) = A2e(A_ideal0);
% stateOrb0(13:18) = ySun;


% Properties
I = computeMOI('res/mass.csv',cm);
[rot,~] = eig(I);

% Run numerical method
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,stateOrb] = ode113(@(t,stateOrb) orbitSimple(t,stateOrb), ...
    tspan,stateOrb0,options);

[t, stateEul] = ode113(@(t,stateEul) kinEulerAngle(t,stateEul,Ix,Iy,Iz), ...
    tspan,stateEul0,options);

% Get Euler angles of ideal rotation
eulerAngs_ideal = nan(size(stateOrb(:,1:3)));
A_ECI2Ideal = nan([3, 3, length(t)]);

for n = 1:length(t)
    r = stateOrb(n,1:3);
    v = stateOrb(n,4:6);
    h = cross(r,v);
    radial = r' / norm(r);
    normal = h' / norm(h);
    tangential = cross(normal, radial);
    A_ECI2RTN = [radial tangential normal]';
    A_ECI2Ideal(:,:,n) = A_RTN2ideal * A_ECI2RTN;
    eulerAngs_ideal(n,:) = A2e(A_ECI2Ideal(:,:,n));
end
eulerAngs_ideal = unwrap(eulerAngs_ideal);

% Get Euler angles of principal axis
eulerAngs_actual = stateEul(:,1:3);
A_ECI2P = nan(3,3,length(t));
for n = 1:length(t)
    A_ECI2P(:,:,n) = e2A(eulerAngs_actual(n,:));
end
eulerAngs_actual = unwrap(eulerAngs_actual);

eulerAngs_error = nan(size(eulerAngs_actual));
for n = 1:length(t)
    A_error = A_RTN2ideal * A_ECI2Ideal(:,:,n);
    eulerAngs_error(n,:) = A2e(A_error);
end

% Plot
figure(1)
hold on
% plot(t/3600, rad2deg(eulerAngs_error))
plot(t/3600, rad2deg(eulerAngs_ideal - eulerAngs_actual))
legend(["\phi_{error}", "\theta_{error}", "\psi_{error}"])
xlabel("time [hr]"); ylabel("Euler Angles [deg]")
% xlim([0, 24])
% saveAsBool(gcf, 'Images/ps6_problem2.png', savePlots)
% saveAsBool(gcf, 'Images/ps6_problem3.png', savePlots)
hold off

figure(2)
plot(t, rad2deg(eulerAngs_ideal))
legend(["\phi", "\theta", "\psi"])

figure(3)
plot(t, rad2deg(eulerAngs_actual))
legend(["\phi", "\theta", "\psi"])
