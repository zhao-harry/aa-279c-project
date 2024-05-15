function [eulerAngs, omega] = kinEulerStepRK4(eulerAngs, omega, Ix, Iy, Iz, tStep)
    state = zeros(6,1);
    state(1:3) = eulerAngs;
    state(4:6) = omega;

    t = 0;

    k1 = kinEulerAngle(t,state,Ix,Iy,Iz);
    k2 = kinEulerAngle(t+tStep/2,state+(k1*tStep/2),Ix,Iy,Iz);
    k3 = kinEulerAngle(t+tStep/2,state+(k2*tStep/2),Ix,Iy,Iz);
    k4 = kinEulerAngle(t+tStep,state+(k3*tStep),Ix,Iy,Iz);
    nextState = state + tStep * (k1/6 + k2/3 + k3/3 + k4/6);

    eulerAngs = nextState(1:3);
    omega = nextState(4:6);
end