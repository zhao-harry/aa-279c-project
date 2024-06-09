close all; clear; clc;
savePlots = false;
modelVars

%% Plotting
time = squeeze(out.q.time);
Lw = squeeze(out.Lw.Data);
w_wheel = Lw ./ control.IWheel;

figure()
for n = 1:4
    subplot(2,2,n)
    hold on
    title(sprintf("Wheel %i", n))
    
    yyaxis left
    ylabel("Angular Velocity [rad/s]")
    plot(time, w_wheel(n,:), 'b')
    ylim([min(w_wheel(:)), max(w_wheel(:))])
    xlabel("Time [s]")
end

%%
syms w0 A11 A12 A13 A21 A22 A23 A31 A32 A33
A = [A11 A12 A13;
        A21 A22 A23;
        A31 A32 A33]
w = [0; -w0; 0];

-crossMatrix(w.'*A.')

%%
wEstimated = w0;
B = [4.869e-6, 7.704e-6, -2.123e-5];
Lw = [0; 0; 0; 0];
Mc = [-0.004307; -0.016; 0.01978];

h = control.A * Lw;
W = [0 0 wEstimated(2);
        0 0 0;
        -wEstimated(2) 0 0];
% W = crossMatrix(W);
Q = -crossMatrix(B);
F = control.F;
R = control.R;
% K = icare(W, Q, F, R)
K = care(W, Q, F, R);
g = inv(W' - K*Q*inv(R)*Q') * K*Mc;

D = - inv(R) * Q' * (K*h - g)

% Limit D if needed
Mmag = cross(D, B);

%% Momentum wheel investigation

Lw_data = out.Lw.data;
LwDot_data = out.LwDot.data;

for n = 100
    Lw = Lw_data(:,:,n);
    Lwdot = Lw_data(:,:,n);
    Lwdot(4) = - abs(Lwdot(4));
    Lw(4) = 51;
    
    saturated = (abs(Lw) > control.LMaxWheel) == 1;
    LwIncreasing = (sign(Lwdot) == sign(Lw));
    Lwdot(saturated & LwIncreasing) = 0;

end