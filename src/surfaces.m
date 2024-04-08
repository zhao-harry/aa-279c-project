function [barycenter, normal, area] = surfaces(filename)
    data = readmatrix(filename);
    barycenter = data(:,1:3);
    normal = data(:,4:6);
    area = data(:,7);
end