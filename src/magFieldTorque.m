function M = magFieldTorque(m, R, t, RE, UT1)
    % Calculated the expected torque due to Earth's magnetic field
    % Inputs:
    % - m: magnetic moment of satellite (N*m/T)
    % - R: position vector of satellite (km)
    % - t: time of simulation (seconds)
    % - RE: radius of Earth (km) (6378 km)
    % - UT1: start time of simulation

    % Find lambda (longitude) and theta (colatitude)
    GMST = time2GMST(t, UT12MJD(UT1));
    [lambda, phi, ~] = ECEF2Geoc(ECI2ECEF(R, GMST));
    theta = pi/2 - phi;

    B = magFieldEarth()
end