function A = euler2A(euler)
    phi = euler(1);
    theta = euler(2);
    psi = euler(3);

    A11 = cos(psi)*cos(phi) - cos(theta)*sin(psi)*sin(phi);
    A12 = cos(psi)*sin(phi) + cos(theta)*sin(psi)*cos(phi);
    A13 = sin(theta)*sin(psi);
    A21 = -sin(psi)*cos(phi) - cos(theta)*cos(psi)*sin(phi);
    A22 = -sin(psi)*sin(phi) + cos(theta)*cos(psi)*cos(phi);
    A23 = sin(theta)*cos(psi);
    A31 = sin(theta)*sin(phi);
    A32 = -sin(theta)*cos(phi);
    A33 = cos(theta);

    A = [A11, A12, A13;
            A21, A22, A23;
            A31, A32, A33];
end