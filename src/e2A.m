function A = e2A(e)
    phi = e(1);
    theta = e(2);
    psi = e(3);
    A11 = cos(phi) * cos(psi) - sin(phi) * cos(theta) * sin(psi);
    A12 = cos(phi) * sin(psi) + sin(phi) * cos(theta) * cos(psi);
    A13 = sin(phi) * sin(theta);
    A21 = -sin(phi) * cos(psi) - cos(phi) * cos(theta) * sin(psi);
    A22 = -sin(phi) * sin(psi) + cos(phi) * cos(theta) * cos(psi);
    A23 = cos(phi) * sin(theta);
    A31 = sin(theta) * sin(psi);
    A32 = -sin(theta) * cos(psi);
    A33 = cos(theta);
    A = [A11 A12 A13; ...
        A21 A22 A23; ...
        A31 A32 A33];
end