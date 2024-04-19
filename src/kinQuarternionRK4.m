function [q,w] = kinQuarternionRK4(q0,w0,Ix,Iy,Iz,tFinal,tStep)
    % 4th order Runge-Kutta integration for quaternions, angular velocity
    nStep = ceil(tFinal/tStep);
    q = nan(nStep+1,4);
    w = nan(nStep+1,3);
    q(1,:) = q0';
    w(1,:) = w0';
    for i = 1:nStep
        t = i * tStep;
        qi = q(i,:)';
        wi = w(i,:)';
        state = [qi;wi];
        k1 = kinQuaternion(t,state,Ix,Iy,Iz);
        k2 = kinQuaternion(t+tStep/2,state+(k1*tStep/2),Ix,Iy,Iz);
        k3 = kinQuaternion(t+tStep/2,state+(k2*tStep/2),Ix,Iy,Iz);
        k4 = kinQuaternion(t+tStep,state+(k3*tStep),Ix,Iy,Iz);
        nextState = state + tStep * (k1/6 + k2/3 + k3/3 + k4/6);
        q(i+1,:) = nextState(1:4) / norm(nextState(1:4));
        w(i+1,:) = nextState(5:7);
    end
end