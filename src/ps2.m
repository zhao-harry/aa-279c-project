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
xlabel('i [km]')
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
Ix = I(1,1);
Iy = I(2,2);
Iz = I(3,3);
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
set(quiver,'XData',[0;0;0])
set(quiver,'YData',[0;0;0])
set(quiver,'ZData',[0;0;0])
set(textx,'Position',[4 0 0])
set(texty,'Position',[0 4 0])
set(textz,'Position',[0 0 4])

quiverPrincipal = copyobj(quiver,gca);
textxPrincipal = copyobj(textx,gca);
textyPrincipal = copyobj(texty,gca);
textzPrincipal = copyobj(textz,gca);
set(quiver,'Color',[0 1 0])
set(quiver,'UData',4.14 * rot(1,:)')
set(quiver,'VData',4.14 * rot(2,:)')
set(quiver,'WData',4.14 * rot(3,:)')
set(quiver,'XData',repmat(cm(1),3,1))
set(quiver,'YData',repmat(cm(2),3,1))
set(quiver,'ZData',repmat(cm(3),3,1))
set(textx,'String','x''')
set(texty,'String','y''')
set(textz,'String','z''')
set(textx,'Position',4 * xPrincipal + cm)
set(texty,'Position',4 * yPrincipal + cm)
set(textz,'Position',4 * zPrincipal + cm)
saveas(gcf,'Images/ps2_model.png');

%% Problem 5
% Define w
w0 = deg2rad([8;4;6]);

tspan = 0:120;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,w] = ode113(@(t,w) eulerPropagator(t,w,Ix,Iy,Iz),tspan,w0,options);

w = rad2deg(w);

figure(1)
plot(t,w)
legend('\omega_{x}','\omega_{y}','\omega_{z}','Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) ']'])
saveas(1,'Images/ps2_euler_equations.png')

%% Problem 6 (Energy)
wx0 = w0(1);
wy0 = w0(2);
wz0 = w0(3);
syms wx wy wz

% Energy ellipsoid
T = (wy0^2*Iy + wz0^2*Iz + wx0^2*Ix) / 2;
energyEllipsoid = wy.^2/(2*T/Iy) + wz.^2/(2*T/Iz) + wx.^2/(2*T/Ix) - 1;
wzEnergy = solve(energyEllipsoid,wz);
fsurf(wzEnergy)
hold on

% Momentum ellipsoid
L = sqrt(wy0^2*Iy^2 + wz0^2*Iz^2 + wx0^2*Ix^2);
momentumEllipsoid = wy.^2/(L^2/Iy^2) + wz.^2/(L^2/Iz^2) + wx.^2/(L^2/Ix^2) - 1;
wzMomentum = solve(momentumEllipsoid,wz);
fsurf(wzMomentum)
hold off