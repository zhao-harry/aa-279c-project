function cm = computeCM(filename)
    data = readmatrix(filename);
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    m = data(:,4);
    cm = [dot(x,m); ...
        dot(y,m); ...
        dot(z,m)] / sum(m);
end