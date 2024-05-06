function M = magFieldTorque(m, R, RE)
    % Calculated the expected torque due to Earth's magnetic field
    % Inputs:
    % - m: magnetic moment of satellite
    % - R: position vector of satellite
    % - RE: radius of Earth

    % Find gamma (longitude) and lambda (colatitude)
    R_ECEF = R;
    R_LLA = functiontoexist(R);
    theta = pi/2 - phi;

    % Find B
    R_norm = norm(R);
    dR = R_norm * 0.01;
    dLambda = lambda * 0.01;
    dTheta = theta * 0.01;
    B = - gradMFPot(R_norm, lambda, theta, RE, dR, dLambda, dTheta);

    M = cross(m, B);
end