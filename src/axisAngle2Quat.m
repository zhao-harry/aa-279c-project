function q = axisAngle2Quat(e, phi)
    % NOTE: phi must be in radians
    q = zeros(4,1);
    q(1) = e(1)*sin(phi/2);
    q(2) = e(2)*sin(phi/2);
    q(3) = e(3)*sin(phi/2);
    q(4) = cos(phi/2);
end