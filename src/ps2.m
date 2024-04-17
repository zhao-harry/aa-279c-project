%% Problem Set 2
clear; close all; clc;

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
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,y] = ode113(@propagator,tspan,yECI,options);

plot3(y(:,1),y(:,2),y(:,3),'LineWidth',2,'Color','green')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal
hold on
[xE,yE,zE] = ellipsoid(0,0,0,6378.1,6378.1,6378.1,20);
surface(xE,yE,zE,'FaceColor','blue','EdgeColor','black');
hold off
saveas(gcf,'Images/ps2_problem1.png');

%% Problem 2
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);
xPrincipal = rot(:,1);
yPrincipal = rot(:,2);
zPrincipal = rot(:,3);

%% Problem 3
figure
gm = importGeometry('res/NISAR.stl');
pdegplot(gm);

quiver = findobj(gca,'type','Quiver');
textx = findobj(gca,'type','Text','String','x');
texty = findobj(gca,'type','Text','String','y');
textz = findobj(gca,'type','Text','String','z');
set(quiver,"XData",[0;0;0])
set(quiver,"YData",[0;0;0])
set(quiver,"ZData",[0;0;0])
set(textx,"Position",[4 0 0])
set(texty,"Position",[0 4 0])
set(textz,"Position",[0 0 4])

quiverPrincipal = copyobj(quiver,gca);
textxPrincipal = copyobj(textx,gca);
textyPrincipal = copyobj(texty,gca);
textzPrincipal = copyobj(textz,gca);
set(quiver,"Color",[0 1 0])
set(quiver,"UData",4.14 * rot(1,:)')
set(quiver,"VData",4.14 * rot(2,:)')
set(quiver,"WData",4.14 * rot(3,:)')
set(quiver,"XData",repmat(cm(1),3,1))
set(quiver,"YData",repmat(cm(2),3,1))
set(quiver,"ZData",repmat(cm(3),3,1))
set(textx,"String",'x''')
set(texty,"String",'y''')
set(textz,"String",'z''')
set(textx,"Position",4 * xPrincipal + cm)
set(texty,"Position",4 * yPrincipal + cm)
set(textz,"Position",4 * zPrincipal + cm)
saveas(gcf,'Images/ps2_model.png');

%% Problem 5
w0Deg = [8;4;6];
w0 = deg2rad(w0Deg);
tspan = 0:120;
w = eulerPropagator(w0,Ix,Iy,Iz,tspan,'Images/ps2_euler_equations.png');

%% Problem 6
[XE,YE,ZE] = ellipsoidEnergy(IPrincipal,w0,'Images/ps2_problem6_energy.png');
[XM,YM,ZM] = ellipsoidMomentum(IPrincipal,w0,'Images/ps2_problem6_momentum.png');

%% Problem 7
w = polhode(XE,YE,ZE,XM,YM,ZM,w,'Images/ps2_problem7.png');

%% Problem 8
w = polhode2D(w,'none','Images/ps2_problem8.png');

%% Problem 9, x-axis
axis = 'x';
w0Deg = 8*[1;0;0];
w0 = deg2rad(w0Deg);
tspan = 0:120;
marker = 'o';

%% Problem 9, y-axis
axis = 'y';
w0Deg = 8*[0.01;1;0.01];
w0 = deg2rad(w0Deg);
tspan = 0:1200;
marker = 'none';

%% Problem 9, z-axis
axis = 'z';
w0Deg = 8*[0.01;0;1];
w0 = deg2rad(w0Deg);
tspan = 0:120;
marker = 'none';

%% Problem 9
w = eulerPropagator(w0,Ix,Iy,Iz,tspan,['Images/ps2_problem9_euler_equations_', axis, '.png']);

[XE,YE,ZE] = ellipsoidEnergy(IPrincipal,w0,['Images/ps2_problem9_energy_', axis, '.png']);
[XM,YM,ZM] = ellipsoidMomentum(IPrincipal,w0,['Images/ps2_problem9_momentum_', axis, '.png']);

w = polhode(XE,YE,ZE,XM,YM,ZM,w,['Images/ps2_problem9_p7_', axis, '.png']);

w = polhode2D(w,marker,['Images/ps2_problem9_p8_', axis, '.png']);