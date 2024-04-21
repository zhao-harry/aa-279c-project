clear; close all; clc

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

%% Non-Axisymmetric Satellite
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);
xPrincipal = rot(:,1);
yPrincipal = rot(:,2);
zPrincipal = rot(:,3);

%% Problem 6 (Quaternions)
axang0 = [sqrt(1/2) sqrt(1/2) 0 pi/4];
q0 = axang2quat(axang0).';
tFinal = 600;
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

tFinal = 600;
tStep = 0.1;
t = 0:tStep:tFinal;

% Forward Euler
% state = kinEulerAngleForwardEuler(state0,Ix,Iy,Iz,tFinal,tStep);

% ode113
tspan = 0:tStep:tFinal;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode113(@(t,state) kinEulerAngle(t,state,Ix,Iy,Iz), ...
    tspan,state0,options);

eulerAngle = wrapTo360(rad2deg(state(:,1:3)));

figure(3)
plot(t,eulerAngle,'LineWidth',1)
legend('\phi','\theta','\psi', ...
    'Location','southwest')
xlabel('Time [s]')
ylabel('Euler Angle [deg]')
saveas(3,'Images/ps3_problem6_euler.png')

%% Problem 7(a)
% Part a: Angular momentum
tLen = length(t);
L_principal = [Ix Iy Iz] .* w;
L_inertial = nan(size(L_principal));
L_norm = nan(1, tLen);

% Part b: Herpolhode
w_inertial = nan(size(w));

for i = 1:tLen
    % Get rotation matrix
    qi = q(i,:);
    A = q2A(qi);

    % Angular momentum
    L_inertial(i,:) = A' * L_principal(i,:)';
    L_norm(i) = norm(L_inertial(i,:));

    % Angular velocity
    w_inertial(i,:) = A' * w(i,:)';
end

figure(4)
hold on
plot(t, L_inertial)
plot(t, L_norm, 'k--')
xlabel('Time [s]')
ylabel('Angular momentum [kg m^{2}/s]')
legend("L_{1}", "L_{2}", "L_{3}", "||L||")
hold off
saveas(4, 'Images/ps3_problem7a.png')

%% Problem 7(b)
figure(5)
plot3(w_inertial(:,1), w_inertial(:,2), w_inertial(:,3), 'r')
grid on
hold on
quiver3(0, 0, 0, ...
    L_inertial(1,1), L_inertial(1,2), L_inertial(1,3), ...
    1e-4)
quiver3(0, 0, 0, ...
    w_inertial(1,1), w_inertial(1,2), w_inertial(1,3), ...
    1)
xlabel('\omega_{x} [rad/s]')
ylabel('\omega_{y} [rad/s]')
zlabel('\omega_{z} [rad/s]')
legend('Herpolhode', ...
    'Angular momentum (L)', ...
    'Angular velocity (\omega)', ...
    'Location','northwest')
hold off
saveas(5, 'Images/ps3_problem7b.png')

%% For fun kinda thing
saveGif = true;
tGif = 240 / tStep;

L_unit = nan(size(L_inertial));
w_unit = nan(size(w_inertial));
if saveGif == true
    gif = figure;
    for i = 1:tLen
        w_unit(i,:) = w_inertial(i,:)./norm(w_inertial(i,:));
        L_unit(i,:) = L_inertial(i,:)./norm(L_inertial(i,:));
    end

    for i = 1:20:tGif
        plot3(w_unit(1:i,1), w_unit(1:i,2), w_unit(1:i,3), 'r')
        grid on
        hold on
        quiver3(0, 0, 0, w_unit(i,1), w_unit(i,2), w_unit(i,3),1)
        quiver3(0, 0, 0, L_unit(i,1), L_unit(i,2), L_unit(i,3),1)
        hold off
        xlim([-1 1])
        ylim([-1 1])
        zlim([-1 1])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Note: all vectors are normalized')
        legend('Herpolhode','\omega','L','Location','northeast')
        exportgraphics(gif,'Images/ps3_problem7b.gif','Append',true);
    end
end

%% Problem 7(c)
% Generate orbit
a = 7125.48662; % km
e = 0.0011650;
i = 98.40508; % degree
O = -19.61601; % degree
w = 89.99764; % degree
nu = -89.99818; % degree

days = 0.069;
tFinal = days * 86400;
tStep = 1;
tspan = 0:tStep:tFinal;

figure(6)
[t,y] = plotECI(a,e,i,O,w,nu,tspan);
hold on
figure(7)
plotECI(a,e,i,O,w,nu,tspan);
hold on
figure(8)
plotECI(a,e,i,O,w,nu,tspan);
hold on

[q,w] = kinQuaternionRK4(q0,w0,Ix,Iy,Iz,tFinal,tStep);

tLen = length(t);
for i = 1:500:tLen
    % Get rotation matrix
    qi = q(i,:);
    A = q2A(qi);
    % Body axes
    B = rot * A * rot';
    % Position
    pos = y(i,1:3);
    radial = pos / norm(pos);
    tangential = y(i,4:6) / norm(y(i,4:6));
    normal = cross(radial,tangential);
    RTN = [radial' tangential' normal'];
    figure(6);
    plotTriad(gca,pos,A,1e3,'r');
    figure(7);
    plotTriad(gca,pos,B,1e3,'m');
    figure(8);
    plotTriad(gca,pos,RTN,1e3,'b');
end
figure(6);
legend('Orbit','Earth','Principal (x-axis)','Location','northwest')
hold off
saveas(gcf,'Images/ps3_problem7c_principal.png');
figure(7);
legend('Orbit','Earth','Body (x-axis)','Location','northwest')
hold off
saveas(gcf,'Images/ps3_problem7c_body.png');
figure(8);
legend('Orbit','Earth','RTN (radial axis)','Location','northwest')
hold off
saveas(gcf,'Images/ps3_problem7c_rtn.png');