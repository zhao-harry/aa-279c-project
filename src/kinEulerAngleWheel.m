function stateDot = kinEulerAngleWheel(t,state,M,r,Ix,Iy,Iz,Ir)
    % Computes state derivative for Euler angles, angular velocity
    % Adds momentum wheel
    % Assign variables
    phi = state(1);
    theta = state(2);
    w = state(4:7);

    stateDot = zeros(7,1);
    % Angular velocity time derivatives
    wDot = eulerEquationWheel(t,w,M,r,Ix,Iy,Iz,Ir);
    stateDot = zeros(7,1);
    stateDot(4) = wDot(1);
    stateDot(5) = wDot(2);
    stateDot(6) = wDot(3);
    stateDot(7) = wDot(4);
    % Euler angle time derivatives
    % 312
    EPrimeInv = [sin(phi)*sin(theta) cos(phi)*sin(theta) cos(theta); ...
        cos(theta)*cos(phi) -sin(phi)*cos(theta) 0; ...
        sin(phi) cos(phi) 0] * (1 / cos(theta));
    % 313
    % EPrimeInv = [-sin(phi)*cos(theta) -cos(phi)*cos(theta) sin(theta); ...
    %     cos(phi)*sin(theta) -sin(phi)*sin(theta) 0; ...
    %     sin(phi) cos(phi) 0] * (1 / sin(theta));
    stateDot(1:3) = EPrimeInv * w(1:3);
end