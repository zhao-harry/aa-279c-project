function eulerAngs = kinEulerAngleStepper(w0, phi, theta, psi, tf, dt, Ix, Iy, Iz)
    % state = [omega_vector; phi; theta; psi];
    
    % Initialize
    wx = w0(1);
    wy = w0(2);
    wz = w0(3);
    numSteps = ceil(tf/dt);
    eulerAngs = nan(3,numSteps+1);
    eulerAngs(:,1) = [phi; theta; psi];
    t = 0;

    for n = 1:numSteps
        % Euler Equation
        wxDot = (Iy - Iz) / Ix * wy * wz;
        wyDot = (Iz - Ix) / Iy * wz * wx;
        wzDot = (Ix - Iy) / Iz * wx * wy;
    
        % Euler Angles
        phiDot = (wx*sin(psi) + wy*cos(psi))/sin(theta);
        thetaDot = wx*cos(psi) - wy*sin(psi);
        psiDot = wz - (wx*sin(psi) + wy*cos(psi))*cot(theta);
        
        % Update and iterate
        wx = wx + wxDot*dt;
        wy = wy + wyDot*dt;
        wz = wz + wzDot*dt;
        phi = phi + phiDot*dt;
        theta = theta + thetaDot*dt;
        psi = psi + psiDot*dt;

        eulerAngs(:,n+1) = [phi; theta; psi];
        t = t + dt;
    end
end