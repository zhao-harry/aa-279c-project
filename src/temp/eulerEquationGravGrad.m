function stateDot = eulerEquationGravGrad(t,state,Ix,Iy,Iz,c,n)
    % Get c vector
    cx = c(1);
    cy = c(2);
    cz = c(3);

    % Computes state derivatives for angular velocity with gravity gradient
    wx = state(1);
    wy = state(2);
    wz = state(3);

    % Angular time velocities
    stateDot(1) = (Iz-Iy)/Ix * (3*n^2*cy*cz - wy*wz);
    stateDot(2) = (Ix-Iz)/Iy * (3*n^2*cz*cx - wz*wx);
    stateDot(3) = (Iy-Ix)/Iz * (3*n^2*cx*cy - wx*wy);
end