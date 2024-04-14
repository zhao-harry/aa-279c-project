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
% Define w
problem9 = false;
LogicalStr = {'false', 'true'};
fprintf("Problem 9 bool is: %s\n", LogicalStr{problem9+1})

if ~problem9
    w0Deg = [8;4;6];
else
    w0Deg = [1;0;0]*8;
end
w0 = deg2rad(w0Deg);

tspan = 0:120;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,w] = ode113(@(t,w) eulerPropagator(t,w,Ix,Iy,Iz),tspan,w0,options);

wDeg = rad2deg(w);

figure(1)
plot(t,wDeg,'LineWidth',2)
legend('\omega_{x}','\omega_{y}','\omega_{z}','Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])

if ~problem9
    prob5Name = 'Images/ps2_euler_equations.png';
else
    prob5Name = 'Images/ps2_problem9_euler_equations.png';
end
saveas(1,prob5Name)

%% Problem 6
T = sum(IPrincipal * w0.^2,"all") / 2;
L = sqrt(sum((w0.*IPrincipal).^2,"all"));

[XE,YE,ZE] = ellipsoid(0,0,0,sqrt(2*T/Ix),sqrt(2*T/Iy),sqrt(2*T/Iz),50);
[XM,YM,ZM] = ellipsoid(0,0,0,L/Ix,L/Iy,L/Iz,50);

energyEllipsoid = surf(XE,YE,ZE,'FaceAlpha',0.5,'FaceColor','blue','DisplayName','Energy Ellipsoid');
xlabel('\omega_{x} [rad/s]')
ylabel('\omega_{y} [rad/s]')
zlabel('\omega_{z} [rad/s]')
axis equal
hold on
momentumEllipsoid = surf(XM,YM,ZM,'FaceAlpha',0.5,'FaceColor','green','DisplayName','Momentum Ellipsoid');
legend('Location','northwest')
hold off

if ~problem9
    prob6Name = 'Images/ps2_problem6.png';
else
    prob6Name = 'Images/ps2_problem9_p6.png';
end
saveas(1,prob6Name)

%% Problem 7
energyEllipsoid = surf(XE,YE,ZE,'FaceAlpha',0.5,'FaceColor','blue','DisplayName','Energy Ellipsoid');
xlabel('\omega_{x} [rad/s]')
ylabel('\omega_{y} [rad/s]')
zlabel('\omega_{z} [rad/s]')
axis equal
hold on
momentumEllipsoid = surf(XM,YM,ZM,'FaceAlpha',0.5,'FaceColor','green','DisplayName','Momentum Ellipsoid');
plot3(w(:,1),w(:,2),w(:,3),'LineWidth',2,'Color','red','DisplayName','Polhode')
legend('Location','northwest')
hold off

if ~problem9
    prob7Name = 'Images/ps2_problem7.png';
else
    prob7Name = 'Images/ps2_problem9_p7.png';
end
saveas(1,prob7Name)

%% Problem 8
subplot(1,3,1)
plot(w(:,2),w(:,3),'Marker','o')
title('Polhode (along x-axis)')
xlabel('\omega_{y} [rad/s]')
ylabel('\omega_{z} [rad/s]')
axis equal

subplot(1,3,2)
plot(w(:,1),w(:,3),'Marker','o')
title('Polhode (along y-axis)')
xlabel('\omega_{x} [rad/s]')
ylabel('\omega_{z} [rad/s]')
axis equal

subplot(1,3,3)
plot(w(:,1),w(:,2),'Marker','o')
title('Polhode (along z-axis)')
xlabel('\omega_{x} [rad/s]')
ylabel('\omega_{y} [rad/s]')
axis equal

if ~problem9
    prob8Name = 'Images/ps2_problem8.png';
else
    prob8Name = 'Images/ps2_problem9_p8.png';
end
saveas(1,prob8Name)
