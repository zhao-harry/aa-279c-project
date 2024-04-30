function e = A2e(A)
    phi = atan2(A(1,3),A(2,3));
    theta = acos(A(3,3));
    psi = atan2(A(3,1),-A(3,2));
    e = [phi;theta;psi];
end