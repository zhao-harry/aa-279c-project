function [w, eulerAngs, M_gg, y] = eulerEquationGravGradRK4(w0,eulerAng0,Ix,Iy,Iz,orbitParams,tFinal,tStep)
    % 4th order Runge-Kutta integration for angular velocity with a
    % momentum wheel

    % Initialize key loop parameters
    nStep = ceil(tFinal/tStep);
    w = nan(nStep+1,3);
    eulerAngs = nan(nStep+1,3);
    y = nan(nStep+1,6);
    M_gg = nan(nStep+1,3);
    eulerAngs(1,:) = eulerAng0';
    w(1,:) = w0';
    
    % Initialize orbit state
    yECI = oe2eci(orbitParams.a,orbitParams.e,orbitParams.i,orbitParams.O,orbitParams.w,orbitParams.nu);
    y(1,:) = yECI;
    mu_E = 3.986e5; %km^3/s^2
    P = 2*pi*sqrt(orbitParams.a^3/mu_E); %s
    n = 2*pi/P; %rad/s
    R_ECI = yECI(1:3);
    A_E2P = e2A(eulerAngs(1,1:3));
    R_principal = A_E2P * R_ECI;
    c = R_principal / norm(R_principal);
    M = gravGradTorque(Ix,Iy,Iz,c,n);
    M_gg(1,:) = M;

    for i = 1:nStep
        t = i * tStep;

        % Update state with Euler Equations
        state = [eulerAngs(i,:)'; w(i,:)'];
        k1 = kinEulerAngleTorque(t,state,M,Ix,Iy,Iz);
        k2 = kinEulerAngleTorque(t+tStep/2,state+(k1*tStep/2),M,Ix,Iy,Iz);
        k3 = kinEulerAngleTorque(t+tStep/2,state+(k2*tStep/2),M,Ix,Iy,Iz);
        k4 = kinEulerAngleTorque(t+tStep,state+(k3*tStep),M,Ix,Iy,Iz);
        nextState = state + tStep * (k1/6 + k2/3 + k3/3 + k4/6);
        eulerAngs(i+1,:) = nextState(1:3);
        w(i+1,:) = nextState(4:6);

        % Update orbit parameters
        state = yECI;
        k1 = propagator(t, state);
        k2 = propagator(t+tStep/2,state+(k1*tStep/2));
        k3 = propagator(t+tStep/2,state+(k2*tStep/2));
        k4 = propagator(t+tStep,state+(k3*tStep));
        nextState = state + tStep * (k1/6 + k2/3 + k3/3 + k4/6);
        yECI = nextState;
        y(i,:) = yECI;

        % Get c vector
        R_ECI = yECI(1:3);
        A_E2P = e2A(eulerAngs(i,1:3)); %double check this
        R_principal = A_E2P * R_ECI;
        c = R_principal / norm(R_principal);

        % Get gravity gradient torque
        M = gravGradTorque(Ix,Iy,Iz,c,n);
        M_gg(i+1,:) = M;
    end
end