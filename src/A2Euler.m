function euler = A2Euler(A)
    theta = acos(A(3,3));
    phi = -atan2(A(3,1), A(3,2));
    psi = atan2(A(1,3), A(2,3));
    euler = [theta; phi; psi];
end