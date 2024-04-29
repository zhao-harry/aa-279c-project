function stateDot = eulerEquationWheel(t,state,Ix,Iy,Iz,Ir)
    % Compute Euler Equations with a Reaction Wheel
    % Get key values
    wx = state(1);
    wy = state(2);
    wz = state(3);
    wr = state(4);
    
    % Initialize
    stateDot = zeros(4,1);

    % Solve
    a = (Iz - Iy)*wz + Ir*wr;
    b = (Ix - Iz)*wz - Ir*wr;
    stateDot(1) = -a/Ix*wy;
    stateDot(2) = -b/Iy*wx;
end