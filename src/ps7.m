%% Problem 5 (Kalman Filter)
tFinal = 6000;
tStep = 1;
tspan = 1:tStep:tFinal;

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
w0 = [0; -n; 0];
x0 = [q0; w0];

q = q0;
euler = A2e(A_Nominal);
P = eye(3);
for i = tspan
    [qkminus,Pkminus] = timeUpdate(tStep,q(:,i),w0,P);
    q(:,i + 1) = qkminus;
    euler(:,i + 1) = A2e(q2A(qkminus));
    P = Pkminus;
end

figure()
plot(tspan,wrapToPi(euler(:,1:end-1)'))
xlabel('Time [h]')
ylabel('Euler Angle (Principal) [rad]')
legend('\phi','\theta','\psi')
