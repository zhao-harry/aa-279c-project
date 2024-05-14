function rECEF = ECI2ECEF(rECI,GMST)
    alpha = GMST;
    R = [cos(alpha) sin(alpha) 0; ...
        -sin(alpha) cos(alpha) 0; ...
        0 0 1];
    rECEF = R * rECI;
end