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
[t_out,y_out] = ode113(@propagator,tspan,yECI,options);
plot3(y_out(:,1),y_out(:,2),y_out(:,3))
xlabel('i [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal
saveas(gcf,'Images/ps2_problem1.png');

%% Problem 2
cm = computeCM("res/mass.csv");
I = computeMOI("res/mass.csv",cm);

[rot,IPrincipal] = eig(I);
xPrincipal = rot(:,1);
yPrincipal = rot(:,2);
zPrincipal = rot(:,3);

%% Problem 3
figure
gm = importGeometry("res/NISAR.stl");
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