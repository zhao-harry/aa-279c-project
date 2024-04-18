function stateDot = kinQuaternion(t,state,w,Ix,Iy,Iz)
    % state = [omega_vector; quat_vector];
    
    % Angular velocity
    wx = state(1);
    wy = state(2);
    wz = state(3);
    stateDot = zeros(7,1);
    stateDot(1) = (Iy - Iz) / Ix * wy * wz;
    stateDot(2) = (Iz - Ix) / Iy * wz * wx;
    stateDot(3) = (Ix - Iy) / Iz * wx * wy;

    % quaternion
    sigma = [0, -wz, -wy, -wx;
                  -wz, 0, wx, wy;
                  wy, -wx, 0, wz;
                  -wx, -wy, -wz, 0];
    
    q = state(4:end);
    qDot = 0.5 * sigma * q;
    stateDot(4:end) = qDot;
end