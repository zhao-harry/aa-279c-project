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

%% Problem 6 (Quaternions)
axang0 = [sqrt(1/2) sqrt(1/2) 0 pi/4];
q0 = axang2quat(axang0).';
tFinal = 120;
tStep = 0.1;
t = 0:tStep:tFinal;

% Forward Euler
% [q,w] = kinQuaternionForwardEuler(q0,w0,Ix,Iy,Iz,tFinal,tStep);

% RK4
[q,w] = kinQuaternionRK4(q0,w0,Ix,Iy,Iz,tFinal,tStep);

figure(2)
hold on
plot(t,q,'LineWidth',1)
legend('q_{1}','q_{2}','q_{3}','q_{4}', ...
    'Location','Southeast')
xlabel('Time [s]')
ylabel('Quaternion')
hold off
saveas(2,'Images/ps3_problem6_quaternions.png')

%% Problem 6 (Euler Angles)
eulerAngle0 = rotm2eul(axang2rotm(axang0))';
state0 = [eulerAngle0;w0];

tFinal = 120;
tStep = 0.1;
t = 0:tStep:tFinal;

% Forward Euler
% state = kinEulerAngleForwardEuler(state0,Ix,Iy,Iz,tFinal,tStep);

% ode113
tspan = 0:0.1:120;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinEulerAngle(t,state,Ix,Iy,Iz), ...
    tspan,state0,options);

eulerAngle = wrapTo360(rad2deg(state(:,1:3)));

figure(3)
plot(t,eulerAngle,'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','Southwest')
xlabel('Time [s]')
ylabel('Euler Angle [deg]')
saveas(3,'Images/ps3_problem6_euler.png')

%% Problem 7
% Part a: Angular momentum
tLen = length(t);
L_sat = [Ix; Iy; Iz] .* wQuats;
L_inertial = nan(size(L_sat));
L_norm = nan(1, tLen);

% Part b: herpolhode
w_inertial = nan(size(wQuats));

for n = 1:tLen
    % Get rotation matrix
    qT = quats(:,n);
    q = qT(1:3);
    q4 = qT(4);
    qx = [0, -q(3), q(2);
            q(3), 0, -q(1);
            -q(2), q(1), 0];
    A = (q4^2 - norm(q))*eye(3) + 2*(q*q') - 2*q4*qx;
    A = quat2rotm(qT');

    % angular momentum
    L_inertial(:,n) = A' * L_sat(:,n);
    L_norm(n) = norm(L_inertial(:,n));

    % angular velocity
    w_inertial(:,n) = A.' * wQuats(:,n);
end

fprintf("FIX THE A MATRIX!\n")

figure(4)
hold on
plot(t, L_inertial)
plot(t, L_norm, 'k--')
legend("L_{1}", "L_{2}", "L_{3}", "||L||")
hold off
saveas(4, 'Images/ps3_problem7a.png')

figure(5)
hold on
plot3(w_inertial)