function [stateDot] = gravGrad(t,state,Ix,Iy,Iz,n)
    % Orbit position and velocity
    r = state(1:3);
    v = state(4:6);

    % Angular velocity
    w = state(7:9);

    % Euler angles
    phi = state(10);
    theta = state(11);

    stateDot = zeros(12,1);
    stateDot(1:3) = v;
    stateDot(4:6) = (-3.986 * 10^5 / norm(r)^2) * r / norm(r); % km/s^2

    h = cross(r,v);
    radial = r / norm(r);
    normal = h / norm(h);
    tangential = cross(normal,radial);
    A_RTN = [radial tangential normal];
    A_ECI2P = e2A(state(10:12));
    A = A_ECI2P * A_RTN';
    c = A(:,1);
    M = gravGradTorque(Ix,Iy,Iz,n,c);
    stateDot(7) = (M(1) - (Iz - Iy) * w(2) * w(3)) / Ix;
    stateDot(8) = (M(2) - (Ix - Iz) * w(3) * w(1)) / Iy;
    stateDot(9) = (M(3) - (Iy - Ix) * w(1) * w(2)) / Iz;

    % Euler angle time derivatives
    % 312
    EPrimeInv = [sin(phi)*sin(theta) cos(phi)*sin(theta) cos(theta); ...
        cos(theta)*cos(phi) -sin(phi)*cos(theta) 0; ...
        sin(phi) cos(phi) 0] * (1 / cos(theta));
    % 313
    % EPrimeInv = [-sin(phi)*cos(theta) -cos(phi)*cos(theta) sin(theta); ...
    %     cos(phi)*sin(theta) -sin(phi)*sin(theta) 0; ...
    %     sin(phi) cos(phi) 0] * (1 / sin(theta));
    stateDot(10:12) = EPrimeInv * w;
end