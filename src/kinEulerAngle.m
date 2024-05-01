function stateDot = kinEulerAngle(t,state,Ix,Iy,Iz)
    % Computes state derivative for Euler angles, angular velocity
    % Assign variables
    phi = state(1);
    theta = state(2);
    w = state(4:6);

    stateDot = zeros(6,1);
    % Angular velocity time derivatives
    stateDot(4) = (Iy - Iz) / Ix * w(2) * w(3);
    stateDot(5) = (Iz - Ix) / Iy * w(3) * w(1);
    stateDot(6) = (Ix - Iy) / Iz * w(1) * w(2);
    % Euler angle time derivatives
    % 312
    EPrimeInv = [sin(phi)*sin(theta) cos(phi)*sin(theta) cos(theta); ...
        cos(theta)*cos(phi) -sin(phi)*cos(theta) 0; ...
        sin(phi) cos(phi) 0] * (1 / cos(theta));
    % 313
    % EPrimeInv = [-sin(phi)*cos(theta) -cos(phi)*cos(theta) sin(theta); ...
    %     cos(phi)*sin(theta) -sin(phi)*sin(theta) 0; ...
    %     sin(phi) cos(phi) 0] * (1 / sin(theta));
    stateDot(1:3) = EPrimeInv * w;
end