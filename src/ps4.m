close all; clear; clc

%% Problem 1a
% Recalculate moments of inertia
cm = computeCM('res/mass.csv');
I = computeMOI('res/mass.csv',cm);

[rot,IPrincipal] = eig(I);
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

% Find Euler kinematics
tFinal = 60;
tStep = 0.01;
t = 0:tStep:tFinal;

w0 = [0; 0; 1];
q0 = [0; 0; 0; 1];

[q,w] = kinQuaternionRK4(q0,w0,Ix,Iy,Iz,tFinal,tStep);

% Plot
figure(1)
plot(t, q)
xlabel('Time [s]')
ylabel('Quaternion')
legend('q_{1}','q_{2}','q_{3}','q_{4}', ...
    'Location','Southeast')
saveas(1,'Images/ps4_problem1a.png')

%% Problem 1b
a = 7125.48662; % km
e = 0.0011650;
i = 98.40508; % degree
O = -19.61601; % degree
w_deg = 89.99764; % degree
nu = -89.99818; % degree

[~,y] = plotECI(a,e,i,O,w_deg,nu,t);
close all

w_RTN = zeros(size(w));

for n = 1:length(t)
    % Get rotation matrix
    qn = q(n,:);
    A = q2A(qn);
    % Body axes
    B = rot * A * rot';
    % Position
    pos = y(n,1:3);
    radial = pos / norm(pos);
    tangential = y(n,4:6) / norm(y(n,4:6));
    normal = cross(radial,tangential);
    RTN = [radial' tangential' normal'];
end

%% Problem 2
% Initial conditions
perturbation = 0.01;
w0x = [1; perturbation; perturbation];
w0y = [perturbation; 1; perturbation];
w0z = [perturbation; perturbation; 1];

w0Mat = {w0x, w0y, w0z};
q0 = [0; 0; 0; 1];

for n = 1:3
    w0 = w0Mat{n};
    
    [q,w] = kinQuaternionRK4(q0,w0,Ix,Iy,Iz,tFinal,tStep);

    figure(2)
    subplot(2,1,1)
    plot(t, w)
    xlabel('Time [s]')
    ylabel('Angular velocity (rad/s)')
    legend('w_{x}','w_{y}','w_{z}', ...
        'Location','Southeast')

    subplot(2,1,2)
    plot(t, q)
    xlabel('Time [s]')
    ylabel('Quaternion')
    legend('q_{1}','q_{2}','q_{3}','q_{4}', ...
        'Location','Southeast')
    saveas(2, ['Images/ps4_problem2a_', sprintf('%i',n), '.png'])
end

