% Satellite orbit initial conditions
a = 7125.48662; % km
e = 0;
i = 98.40508; % degree
O = -19.61601; % degree
w = 89.99764; % degree
nu = -89.99818; % degree
muE = 3.986e5; % km^3 / s^2
n = sqrt(muE / a^3);

% Compute initial position and attitude
y = oe2eci(a,e,i,O,w,nu);
r0 = y(1:3);
v0 = y(4:6);
h = cross(r0,v0);
radial = r0 / norm(r0);
normal = h / norm(h);
tangential = cross(normal,radial);
A_Nominal = [-radial -normal -tangential]';

q0 = A2q(A_Nominal);
w0 = [-n; -n; -n];
x0 = [q0; w0];

qkminus1 = q0;
wkminus1 = w0;
euler = A2e(A_Nominal);
Pkminus1 = eye(6);
Q = Pkminus1 / 100;

dt = 0.01;
u = [0 0 0]';

yk = [[1:10; 1:10; 1:10] w0];
% yk = [yk, zeros(size(w0)); zeros(size(yk)), w0];
hk = [[1:10; 1:10; 1:10] w0];
% hk = [hk, zeros(size(w0)); zeros(size(hk)), w0];

[qkminus, xkminus, Pkminus] = timeUpdate(dt, qkminus1, wkminus1, Pkminus1, u, constants, sensors);

[qkplus,wkplus,Pkplus] = measurementUpdate(qkminus, Pkminus, xkminus, yk, hk, u, sensors);

%% Function
function [qkminus, xkminus, Pkminus] = timeUpdate(dt, qkminus1, wkminus1, Pkminus1, u, constants, sensors)
    % MEKF, only time update implemented

    % Get values
    Ix = constants.Ix; Iy = constants.Iy; Iz = constants.Iz;
    Q = sensors.Q;
    xkminus1 = [zeros([3 1]); wkminus1];

    w = norm(wkminus1);
    wcross = [0, -wkminus1(3), wkminus1(2); ...
        wkminus1(3), 0, -wkminus1(1); ...
        -wkminus1(2), wkminus1(1), 0];
    cw = cos(w * dt / 2);
    sw = sin(w * dt / 2);

    A = eye(3) + sw * wcross / w + (1 - cw) * wcross^2 / w^2;
    phi = eye(3) + 0.5 * dt * ...
        [0, (Iy - Iz) / Ix * wkminus1(3), (Iy - Iz) / Ix * wkminus1(2); ...
        (Iz - Ix) / Iy * wkminus1(3), 0, (Iz - Ix) / Iy * wkminus1(1); ...
        (Ix - Iy) / Iz * wkminus1(2), (Ix - Iy) / Iz * wkminus1(1), 0];
    B = dt * [1 / Ix, 0, 0; ...
        0, 1 / Iy, 0; ...
        0, 0, 1 / Iz];
    STM = [A, zeros(3); zeros(3), phi];
    O = eye(4) + 0.5 * dt * ...
        [0, wkminus1(3), -wkminus1(2), wkminus1(1); ...
        -wkminus1(3), 0, wkminus1(1), wkminus1(2); ...
        wkminus1(2), -wkminus1(1), 0, wkminus1(3); ...
        -wkminus1(1), -wkminus1(2), -wkminus1(3), 0];

    xkminus = STM * xkminus1 + [zeros([3 1]); B * u];
    wkminus = xkminus(4:6);
    qkminus = O * qkminus1;
    qkminus = qkminus / norm(qkminus);
    Pkminus = STM * Pkminus1 * STM' + Q;
end

function [qkplus,wkplus,Pkplus] = measurementUpdate(qkminus, Pkminus, xkminus, yk, hk, u, sensors)
    % Get R
    numVectMeasurements = size(hk,2)-1;
    R = diag([repelem(sensors.tracker_error, 3*(numVectMeasurements-1)), repelem(sensors.sun_error, 3), repelem(sensors.gyro_error, 3)]).^2;

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
    Kk = Pkminus * Hk' * inv(Hk * Pkminus * Hk' + R);
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