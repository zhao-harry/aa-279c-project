function stateDot = kinWheel(t,state,M,r,Ix,Iy,Iz,Ir)
    % Computes state derivative for Euler angles, angular velocity
    % Adds momentum wheel
    % Assign variables
    theta = state(2);
    psi = state(3);
    wx = state(4);
    wy = state(5);
    wz = state(6);
    wr = state(7);
    w = [wx; wy; wz; wr];

    stateDot = zeros(7,1);
    % Angular velocity time derivatives
    wDot = eulerEquationWheel(t,w,M,r,Ix,Iy,Iz,Ir);
    stateDot = zeros(7,1);
    stateDot(4) = wDot(1);
    stateDot(5) = wDot(2);
    stateDot(6) = wDot(3);
    stateDot(7) = wDot(4);
    % Euler angle time derivatives
    stateDot(1) = (wx*sin(psi) + wy*cos(psi))/sin(theta);
    stateDot(2) = wx*cos(psi) - wy*sin(psi);
    stateDot(3) = wz - (wx*sin(psi) + wy*cos(psi))*cot(theta);
end