function [wDot] = eulerPropagator(t, w, Ix, Iy, Iz)
    wx = w(1);
    wy = w(2);
    wz = w(3);
    wDot(1) = (Iy-Iz)/Ix*wy*wz;
    wDot(2) = (Iz-Ix)/Iy*wz*wx;
    wDot(3) = (Ix-Iy)/Iz*wx*wy;
end