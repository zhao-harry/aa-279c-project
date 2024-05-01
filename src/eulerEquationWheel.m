function [wDot] = eulerEquationWheel(t,w,M,r,Ix,Iy,Iz,Ir)
    wx = w(1); wy = w(2); wz = w(3); wr = w(4);
    Mx = M(1); My = M(2); Mz = M(3); Mr = M(4);
    rx = r(1); ry = r(2); rz = r(3);
    wrDot = Mr / Ir;
    wDot = zeros(4,1);
    wDot(1) = (Mx - (Ir*wrDot*rx) + (Iy - Iz) * wy * wz + Ir * wr * ...
        (wz * ry - wy * rz)) / Ix;
    wDot(2) = (My - (Ir*wrDot*ry) + (Iz - Ix) * wz * wx + Ir * wr * ...
        (wx * rz - wz * rx)) / Iy;
    wDot(3) = (Mz - (Ir*wrDot*rz) + (Ix - Iy) * wx * wy + Ir * wr * ...
        (wy * rx - wx * ry)) / Iz;
    wDot(4) = wrDot;
end