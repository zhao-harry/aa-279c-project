function [statedot] = propagator(t, state)
    r = state(1:3);
    v = state(4:6);
    statedot = zeros(6,1);
    statedot(1:3) = v; % km/s
    statedot(4:6) = (-3.986 * 10^5 / norm(r)^2) * r / norm(r); % km/s^2
end