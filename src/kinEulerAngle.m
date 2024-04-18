function stateDot = kinEulerAngle(t,state,Ix,Iy,Iz)
    % state = [omega_vector; phi; theta; psi];
    
    % Angular velocity
    wx = state(1);
    wy = state(2);
    wz = state(3);
    stateDot = zeros(6,1);
    stateDot(1) = (Iy - Iz) / Ix * wy * wz;
    stateDot(2) = (Iz - Ix) / Iy * wz * wx;
    stateDot(3) = (Ix - Iy) / Iz * wx * wy;

    % Euler Angles
    phi = state(4);
    theta = state(5);
    psi = state(6);
    phiDot = (wx*sin(psi) + wy*cos(psi))/sin(theta);
    thetaDot = wx*cos(psi) - wy*sin(psi);
    psiDot = wz - (wx*sin(psi) + wy*cos(psi))*cot(theta);
    stateDot(4:end) = [phiDot; thetaDot; psiDot];
end