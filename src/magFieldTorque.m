function M = magFieldTorque(m, R, t, RE)
    % Calculated the expected torque due to Earth's magnetic field
    % Inputs:
    % - m: magnetic moment of satellite
    % - R: position vector of satellite
    % - t: time of simulation (seconds)
    % - RE: radius of Earth

    % Find lambda (longitude) and theta (colatitude)
    UT1 = [2025, 1, 1];
    GMST = time2GMST(t, UT12MJD(UT1));
    [lambda, phi, ~] = ECEF2Geoc(ECI2ECEF(R, GMST));
    theta = pi/2 - phi;

    % Find B
    R_norm = norm(R);
    dR = R_norm * 0.01;
    dLambda = lambda * 0.01;
    dTheta = theta * 0.01;
    B = - gradMFPot(R_norm, lambda, theta, RE, dR, dLambda, dTheta);

    M = cross(m, B);
end