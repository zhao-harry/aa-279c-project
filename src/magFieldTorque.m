function [M, B_ECEF] = magFieldTorque(m,R,eulerAngle,t,RE,UT1)
    % Calculated the expected torque due to Earth's magnetic field
    % Inputs:
    % - m: magnetic moment of satellite [N * m / T]
    % - R: position vector of satellite [km]
    % - t: time of simulation [s]
    % - RE: radius of Earth [km] (6378 km)
    % - UT1: start time of simulation

    % Find lambda (longitude) and theta (colatitude)
    GMST = time2GMST(t,UT12MJD(UT1));
    [lat,lon,~] = ECEF2Geoc(ECI2ECEF(R,GMST),t);
    lambda = lat;
    phi = lon;
    theta = pi/2 - phi;

    [B_R,B_theta,B_phi] = magFieldEarth(R,lambda,theta,RE);
%     B_latlon = [B_R, B_theta, B_phi];
    
    timeVec = UT1;
    timeVec(2) = timeVec(2) + t/86400;
    [XYZ,H,D,I,F] = wrldmagm(norm(R)-RE,lat,lon,decyear(timeVec));

    delta = phi;
    alpha = lambda + GMST;

    B_x = (B_R * cos(delta) + B_theta * sin(delta)) * ...
        cos(alpha) - B_phi * sin(alpha);
    B_y = (B_R * cos(delta) + B_theta * sin(delta)) * ...
        sin(alpha) + B_phi * cos(alpha);
    B_z = (B_R * sin(delta) - B_theta * cos(delta));

    B_ECI = [B_x;B_y;B_z];
    B_ECEF = ECI2ECEF(B_ECI, GMST);
    B = e2A(eulerAngle) * B_ECI;

    M = cross(m,B);
end