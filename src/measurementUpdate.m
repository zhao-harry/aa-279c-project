function [qkplus,wkplus,Pkplus,xkplus,z] = ...
    measurementUpdate(qkminus, Pkminus, xkminus, yk, hk, u, sensors)
    % Get R
    numVectMeasurements = size(hk,2)-1;
    R = diag([repelem(sensors.trackerError, 3*(numVectMeasurements-1)), ...
        repelem(sensors.sunError, 3), repelem(sensors.gyroError, 3)]).^2;

    % Compute residual
    z = yk - hk;

    % Adjust yk and hk vectors
    yk = yk(:);
    hk = hk(:);
    
    % Get Hk
    Hk = zeros([3*(numVectMeasurements+1), 6]);
    for n = 1:3:(3*numVectMeasurements)
        Hk_vect = crossMatrix(hk(n:n+2));
        Hk(n:n+2,:) = [Hk_vect, zeros(3)];
    end
    Hk(end-2:end,:) = [zeros(3), eye(3)];
    Kk = (Pkminus * Hk') / (Hk * Pkminus * Hk' + R);
    xkplus = xkminus + Kk * (yk - hk);
    Pkplus = Pkminus - Kk * Hk * Pkminus;
    wkplus = xkplus(4:6);
    ax = xkplus(1); ay = xkplus(2); az = xkplus(3);

    % Attitude update
    qkplus = [1, az/2, -ay/2, ax/2; ...
        -az/2, 1, ax/2, ay/2; ...
        ay/2, -ax/2, 1, az/2; ...
        -ax/2, -ay/2, -az/2, 1] * qkminus;
end