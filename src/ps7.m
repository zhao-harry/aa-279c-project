%% Import mass properties
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

%% Problem 5 (Kalman Filter)
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
w0 = [0.001; -n; 0.001];
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
plot(tspan,euler(:,1:end-1)')
xlabel('Time [s]')
ylabel('Euler Angle (Principal) [rad]')
legend('\phi','\theta','\psi')

figure()
plot(tspan,w(:,1:end-1)')
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
legend('w_{x}','w_{y}','w_{z}')

% Simulation state
state0 = [A2e(A_Nominal); w0];
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinEulerAngle(t,state,Ix,Iy,Iz), ...
    tspan,state0,options);

% Errors
figure()
plot(t,euler(:,1:end-1)' - state(:,1:3))
ylim([-1e-3 1e-3])
xlabel('Time [s]')
ylabel('Euler Angle Error [rad]')
legend('\phi','\theta','\psi')

figure()
plot(t,w(:,1:end-1)' - state(:,4:6))
xlabel('Time [s]')
ylabel('Angular Velocity Error [rad/s]')
legend('w_{x}','w_{y}','w_{z}')