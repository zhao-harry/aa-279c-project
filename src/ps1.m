%% Center of mass
cm = computeCM("res/mass.csv");

%% Surface properties
[barycenter, normal, area] = surfaces("res/area.csv");

%% Plot spacecraft with body axes
figure
gm = importGeometry("res/NISAR.stl");
pdegplot(gm);
saveas(gcf,'Images/ps1_model.png');