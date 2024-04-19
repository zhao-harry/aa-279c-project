function [state] = kinEulerAngleForwardEuler(state0,Ix,Iy,Iz,tFinal,tStep)
    % Forward Euler integration for state Euler angles, angular velocity
    nStep = ceil(tFinal/tStep);
    state = nan(nStep+1,6);
    state(1,:) = state0;
    for i = 1:nStep
        t = i * tStep;
        statei = state(i,:)';
        stateDot = kinEulerAngle(t,statei,Ix,Iy,Iz);
        state(i+1,:) = statei + tStep * stateDot;
    end
end