function w = eulerEquationWheelRK4(w0,Ix,Iy,Iz,tFinal,tStep)
    % 4th order Runge-Kutta integration for angular velocity with a
    % momentum wheel
    nStep = ceil(tFinal/tStep);
    w = nan(nStep+1,3);
    w(1,:) = w0';

    for i = 1:nStep
        t = i * tStep;
        wi = w(i,:)';
        state = wi;
        k1 = eulerEquationWheel(t,state,Ix,Iy,Iz);
        k2 = eulerEquationWheel(t+tStep/2,state+(k1*tStep/2),Ix,Iy,Iz);
        k3 = eulerEquationWheel(t+tStep/2,state+(k2*tStep/2),Ix,Iy,Iz);
        k4 = eulerEquationWheel(t+tStep,state+(k3*tStep),Ix,Iy,Iz);
        nextState = state + tStep * (k1/6 + k2/3 + k3/3 + k4/6);
        w(i+1,:) = nextState;
    end
end