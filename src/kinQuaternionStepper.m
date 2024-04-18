function quats = kinQuaternionStepper(w0, qT, tf, dt, Ix, Iy, Iz)
    % state = [omega_vector; phi; theta; psi];
    
    % Initialize
    wx = w0(1);
    wy = w0(2);
    wz = w0(3);
    numSteps = ceil(tf/dt);
    quats = nan(4,numSteps+1);
    quats(:,1) = qT;
    t = 0;

    for n = 1:numSteps
        % Euler Equation
        wxDot = (Iy - Iz) / Ix * wy * wz;
        wyDot = (Iz - Ix) / Iy * wz * wx;
        wzDot = (Ix - Iy) / Iz * wx * wy;
    
        % Quaternions
        sigma = [0, -wz, -wy, -wx;
                      -wz, 0, wx, wy;
                      wy, -wx, 0, wz;
                      -wx, -wy, -wz, 0];
        qDot = 0.5 * sigma * qT;
        
        % Update and iterate
        wx = wx + wxDot*dt;
        wy = wy + wyDot*dt;
        wz = wz + wzDot*dt;
        qT = qT + qDot*dt;

        quats(:,n+1) = qT;
        t = t + dt;
    end
end