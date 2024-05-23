%%
close all; clear; clc
savePlots = true;

%% Import mass properties
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

%% Problem 1
% NOTE: Requires fixed timestep (0.1s) to get good data
qVals = squeeze(out.q.data);
qMeasVals = squeeze(out.qMeasured.data);
timeVals = squeeze(out.q.time);
timeValsMeas = squeeze(out.qMeasured.time);

eulerVals = quats2Euler(qVals);
eulerValsMeas = quats2Euler(qMeasVals);

figure(1)
hold on
plot(timeValsMeas, rad2deg(eulerValsMeas), 'b')
plot(timeVals, rad2deg(eulerVals), 'r--')
hold off
%%
eulerError = zeros(size(eulerVals));
A_error = zeros(3,3,length(timeVals));
for n = 1:length(timeVals)
    A_ECI2true = e2A(eulerVals(:,n));
    A_ECI2err = e2A(eulerValsMeas(:,n));
    A_true2err = A_ECI2err * A_ECI2true';
    eulerError(:,n) = A2e(A_true2err);
    A_error(:,:,n) = A_true2err;
end

figure(1)
plot(timeVals/3600, rad2deg(eulerError))
xlim([0 timeVals(end)/3600])
xlabel("time [hr]")
ylabel("euler angle [deg]")
legend("\phi", "\theta", "\psi")
ylim([-20, 20])
saveas(1, "Images/ps7_problem2_qMethod.png")
% saveas(1, "Images/ps7_problem2_DAD.png")
% saveas(1, "Images/ps7_problem2_DADFict.png")
% saveas(1, "Images/ps7_problem2_kin.png")

%% Problem 3
frob_norm = nan(1, length(timeVals));

for n = 1:length(timeVals)
    phi = eulerError(1,n);
    theta = eulerError(2,n);
    psi = eulerError(3,n);
    A_smallAng = [1 phi -psi;
                          -phi 1 theta;
                          psi -theta 1];
    frob_norm(n) = norm(A_error(:,:,n) - A_smallAng);
end

figure(2)
plot(timeVals/3600, frob_norm)
xlabel("time [hr]")
ylabel("||A_{error} - A_{small angle}||_F")
xlim([0 timeVals(end)/3600])

% saveas(2, "Images/ps7_problem3.png")

%% Problem 5-6
tFinal = 18000;
tStep = 0.1;
tspan = 1:tStep:tFinal;

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

q0 = A2q(A_Nominal);
w0 = [-n; -n; -n];
x0 = [q0; w0];

q = q0;
w = w0;
euler = A2e(A_Nominal);
P = eye(6);
Q = P / 100;
for i = 1:length(tspan)
    u = zeros([3 1]);
    [q(:,i + 1),w(:,i + 1),P] = timeUpdate(tStep,q(:,i),w(:,i),P,Q,Ix,Iy,Iz,u);
    euler(:,i + 1) = A2e(q2A(q(:,i + 1)));
end

figure()
plot(tspan,wrapToPi(euler(:,1:end-1)'))
xlabel('Time [s]')
ylabel('Euler Angle (Principal) [rad]')
legend('\phi','\theta','\psi')
saveAsBool(gcf,'Images/ps7_problem5a_angle_est.png',savePlots)

figure()
plot(tspan,wrapToPi(w(:,1:end-1)'))
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
legend('w_{x}','w_{y}','w_{z}')
saveAsBool(gcf,'Images/ps7_problem5a_angvel_est.png',savePlots)

% Simulation state
state0 = [A2e(A_Nominal); w0];
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinEulerAngle(t,state,Ix,Iy,Iz), ...
    tspan,state0,options);

figure()
plot(tspan,wrapToPi(state(:,1:3)))
xlabel('Time [s]')
ylabel('Euler Angle (Principal) [rad]')
legend('\phi','\theta','\psi')
saveAsBool(gcf,'Images/ps7_problem5a_angle_sim.png',savePlots)

figure()
plot(tspan,wrapToPi(state(:,4:6)))
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
legend('w_{x}','w_{y}','w_{z}')
saveAsBool(gcf,'Images/ps7_problem5a_angvel_sim.png',savePlots)

% Errors
eulerError = size(state(:,1:3));
for i = 1:length(t)
    A_MEKF = e2A(euler(:,i));
    A_Sim = e2A(state(i,1:3));
    eulerError(i,1:3) = A2e(A_MEKF * A_Sim');
end

figure()
plot(t,wrapToPi(eulerError))
ylim([-0.01 0.01])
xlabel('Time [s]')
ylabel('Euler Angle Error [rad]')
legend('\phi','\theta','\psi')
saveAsBool(gcf,'Images/ps7_problem6_angle_err.png',savePlots)

figure()
plot(t,w(:,1:end-1)' - state(:,4:6))
xlabel('Time [s]')
ylabel('Angular Velocity Error [rad/s]')
legend('w_{x}','w_{y}','w_{z}')
saveAsBool(gcf,'Images/ps7_problem6_angvel_err.png',savePlots)