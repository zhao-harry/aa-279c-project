%% Center of mass
cm = computeCM("res/mass.csv");

%% Moment of inertia
origin = [0;0;0];
I = computeMOI("res/mass.csv",origin);

%% Surface properties
[barycenter,normal,area] = surfaces("res/area.csv");

%% Plot spacecraft with body axes
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
saveas(gcf,'Images/ps1_model.png');