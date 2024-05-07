function [stateDot] = orbitTorque(t,state,Ix,Iy,Iz,CD,Cd,Cs,P,barycenter,normal,area,cm,n)
    warning('off','aero:atmosnrlmsise00:setf107af107aph')

    % Orbit position and velocity
    r = state(1:3);
    v = state(4:6);
    rEarth = state(13:15);
    vEarth = state(16:18);

    % Angular velocity
    w = state(7:9);

    % Euler angles
    phi = state(10);
    theta = state(11);

    % Gravity gradient torque
    radial = r / norm(r);
    A_ECI2P = e2A(state(10:12));
    c = A_ECI2P * radial;
    Mgg = gravGradTorque(Ix,Iy,Iz,n,c);

    % Drag
    [~,density] = atmosnrlmsise00(1000 * (norm(r) - 6378.1),0,0,2000,1,0);
    rho = density(6);
    vPrincipal = A_ECI2P * v;
    [~,Md] = drag(vPrincipal,rho,CD,barycenter,normal,area,cm);

    % Solar radiation pressure
    Rx = [1 0 0; 0 cosd(23.5) -sind(23.5); 0 sind(23.5) cosd(23.5)];
    s = A_ECI2P * (-Rx * rEarth - r); % Double check this
    [~,Msrp] = srp(s,P,Cd,Cs,barycenter,normal,area,cm);

    % Compute net moments
    Mx = Mgg(1) + Md(1) + Msrp(1);
    My = Mgg(2) + Md(2) + Msrp(2);
    Mz = Mgg(3) + Md(3) + Msrp(3);

    % Time derivatives
    stateDot = zeros(12,1);
    stateDot(1:3) = v;
    stateDot(4:6) = (-3.986E5 / norm(r)^2) * r / norm(r); % km/s^2
    stateDot(7) = (Mx - (Iz - Iy) * w(2) * w(3)) / Ix;
    stateDot(8) = (My - (Ix - Iz) * w(3) * w(1)) / Iy;
    stateDot(9) = (Mz - (Iy - Ix) * w(1) * w(2)) / Iz;

    % 312 Euler angle time derivatives
    EPrimeInv = [sin(phi)*sin(theta) cos(phi)*sin(theta) cos(theta); ...
        cos(theta)*cos(phi) -sin(phi)*cos(theta) 0; ...
        sin(phi) cos(phi) 0] * (1 / cos(theta));
    stateDot(10:12) = EPrimeInv * w;

    % Sun position
    stateDot(13:15) = vEarth;
    stateDot(16:18) = (-1.327E11 / norm(rEarth)^2) * ...
        rEarth / norm(rEarth); % km/s^2
end