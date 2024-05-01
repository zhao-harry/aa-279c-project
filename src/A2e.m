function e = A2e(A)
    % 312
    phi = atan2(A(1,2),A(2,2));
    theta = -asin(A(3,2));
    psi = atan2(A(3,1),A(3,3));
    % 313
    % phi = atan2(A(1,3),A(2,3));
    % theta = acos(A(3,3));
    % psi = atan2(A(3,1),-A(3,2));
    e = [phi;theta;psi];
end