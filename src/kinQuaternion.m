function stateDot = kinQuaternion(t,state,Ix,Iy,Iz)
    q = state(1:4);
    wx = state(5);
    wy = state(6);
    wz = state(7);
    stateDot = zeros(7,1);
    stateDot(5) = (Iy - Iz) / Ix * wy * wz;
    stateDot(6) = (Iz - Ix) / Iy * wz * wx;
    stateDot(7) = (Ix - Iy) / Iz * wx * wy;
    sigma = [0, -wz, -wy, -wx; ...
        -wz, 0, wx, wy; ...
        wy, -wx, 0, wz; ...
        -wx, -wy, -wz, 0];
    qDot = 0.5 * sigma * q;
    stateDot(1:4) = qDot;
end