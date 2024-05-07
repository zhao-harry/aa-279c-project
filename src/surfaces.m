function [barycenter, normal, area] = surfaces(filename,rotm)
    data = readmatrix(filename);
    barycenter = rotm * data(:,1:3)';
    normal = rotm * data(:,4:6)';
    area = data(:,7)';
end