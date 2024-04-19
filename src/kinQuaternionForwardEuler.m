function [q,w] = kinQuaternionForwardEuler(q0,w0,Ix,Iy,Iz,tFinal,tStep)
    % Forward Euler integration for quaternions, angular velocity
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
        stateDot = kinQuaternion(t,state,Ix,Iy,Iz);
        nextState = state + tStep * stateDot;
        q(i+1,:) = nextState(1:4) / norm(nextState(1:4));
        w(i+1,:) = nextState(5:7);
    end
end