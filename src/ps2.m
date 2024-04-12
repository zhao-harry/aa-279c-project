%% Problem Set 2
clear; close all; clc

%% Problem 1
a = 7125.48662; % km
e = 0.0011650;
i = 98.40508; % degree
O = -19.61601; % degree
w = 89.99764; % degree
nu = -89.99818; % degree

yECI = oe2eci(a,e,i,O,w,nu);

days = 0.5;
tspan = 0:days*86400;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, y_out] = ode113(@propagator, tspan, yECI, options);
plot3(y_out(:,1),y_out(:,2),y_out(:,3))
xlabel('i [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal
saveas(gcf,'Images/ps2_problem1.png');

%% Problem 2
I = computeMOI("res/mass.csv");

[rotation_matrix, I_principal] = eig(I);
x_principal = rotation_matrix(:,1);
y_principal = rotation_matrix(:,3);
z_principal = rotation_matrix(:,2);

