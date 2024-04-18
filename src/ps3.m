%% Problem 1
IPrincipal = [7707.07451493673 0 0; ...
    0 7707.0745149367 0; ...
    0 0 18050.0227594212];
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

w0Deg = [8;4;6];
w0 = deg2rad(w0Deg);
tspan = 0:0.1:120;
w = eulerPropagator(w0,Ix,Iy,Iz,tspan,'Images/ps3_problem1.png');

%% Problem 2
lambda = w0(3) * (Iz - Iy) / Ix;
wxy = (w0(1) + w0(2) * 1j) * exp(1j * lambda * tspan);
wx = real(wxy);
wy = imag(wxy);
wz = w0(3) * ones(size(wxy));
wAnalytical = [wx',wy',wz'];
wDegAnalytical = rad2deg(wAnalytical);

figure(1)
plot(tspan,wDegAnalytical,'LineWidth',2)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(1,'Images/ps3_problem2.png')

%% Problem 3
error = w - wAnalytical;
plot(tspan,error,'LineWidth',2)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel('Angular velocity (\omega) [rad/s]')
saveas(gcf,'Images/ps3_problem3.png')

%% Problem 6
% Quaternions
% Initial angles
e0 = [1; 0; 0];
phi0 = pi/2;
q0 = axisAngle2Quat(e0, phi0);
% state = [w0; q0];

tf = 60;
dt = 0.0001;
quats = kinQuaternionStepper(w0, q0, tf, dt, Ix, Iy, Iz);

% ode113 section that doesn't work
% % Numerically propogate quaternions
% tspan = 0:60;
% options = odeset('RelTol',1e-6,'AbsTol',1e-9);
% [t,state] = ode113(@(t,w) kinQuaternion(t,state,Ix,Iy,Iz),tspan,state,options);

% quats = state(:,4:end);
figure(2)
hold on
plot(t, quats,'LineWidth',1)
legend('q_{1}','q_{2}','q_{3}', 'q_{4}', ...
    'Location', 'Southwest')
xlabel('Time [s]')
ylabel('Quaternion')
hold off
% saveas(1,filename)

%% 
% Euler Angles
phi0 = deg2rad(5);
theta0 = deg2rad(5);
psi0 = deg2rad(5);
state = [w0; phi0; theta0; psi];

% NOTE: ode113 method doesn't seem to work
% tspan = 0:60;
% options = odeset('RelTol',1e-6,'AbsTol',1e-9);
% [t,state] = ode113(@(t,w) kinEulerAngle(t,state,Ix,Iy,Iz),tspan,state,options);
% 
% eulerAngs = wrapTo360(rad2deg(state(:,4:end)));

% Stepper
tf = 60;
dt = 0.0001;
eulerAngs = kinEulerAngleStepper(w0, phi0, theta0, psi0, tf, dt, Ix, Iy, Iz);
t = 0:dt:tf;
% saveas(1,'Images/ps4_problem6.png')

figure(3)
plot(t, eulerAngs,'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location', 'Southwest')
xlabel('Time [s]')
ylabel('Euler Angle [deg]')
