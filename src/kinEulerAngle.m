function stateDot = kinEulerAngle(t,state,Ix,Iy,Iz)
    % Computes state derivative for Euler angles, angular velocity
    % Assign variables
    theta = state(2);
    psi = state(3);
    wx = state(4);
    wy = state(5);
    wz = state(6);

    stateDot = zeros(6,1);
    % Angular velocity time derivatives
    stateDot(4) = (Iy - Iz) / Ix * wy * wz;
    stateDot(5) = (Iz - Ix) / Iy * wz * wx;
    stateDot(6) = (Ix - Iy) / Iz * wx * wy;
    % Euler angle time derivatives
    stateDot(1) = (wx*sin(psi) + wy*cos(psi))/sin(theta);
    stateDot(2) = wx*cos(psi) - wy*sin(psi);
    stateDot(3) = wz - (wx*sin(psi) + wy*cos(psi))*cot(theta);
end