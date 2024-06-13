function B = getMagField(rECI,A_ECI2P,t,constants)
    % Calculated the expected torque due to Earth's magnetic field
    % Inputs:
    % - m: magnetic moment of satellite [N * m / T]
    % - R: position vector of satellite [km]
    % - t: time of simulation [s]
    % - RE: radius of Earth [km] (6378 km)
    % - UT1: start time of simulation

    GMST = time2GMST(t,UT12MJD(constants.UT1));
    [lat,lon,~] = ECEF2Geoc(ECI2ECEF(rECI,GMST),t);
    theta = pi/2 - lat;

    [B_R,B_theta,B_phi] = magFieldEarth(rECI,lon,theta,constants.RE);

    delta = lat;
    alpha = lon + GMST;

    B_x = (B_R * cos(delta) + B_theta * sin(delta)) * ...
        cos(alpha) - B_phi * sin(alpha);
    B_y = (B_R * cos(delta) + B_theta * sin(delta)) * ...
        sin(alpha) + B_phi * cos(alpha);
    B_z = (B_R * sin(delta) - B_theta * cos(delta));

    B_ECI = [B_x;B_y;B_z];
    B = A_ECI2P * B_ECI;
end