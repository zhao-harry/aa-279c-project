close all; clear; clc;
savePlots = true;
modelVars

%%
q = squeeze(out.q.data);
qMeas = squeeze(out.qMeasured.data);
time = squeeze(out.q.time);
timeMeas = squeeze(out.qMeasured.time);

eulerVals = quats2Euler(q);
eulerValsMeas = quats2Euler(qMeas);

figure()
hold on
plot(time, rad2deg(eulerVals))
legend(["\phi" "\theta" "\psi"])
hold off

%% Problem 2
time = squeeze(out.q.time);
sinWave = squeeze(out.sinWave.Data);
Lw = squeeze(out.Lw.Data);
w_wheel = Lw ./ control.IWheel;

figure()
for n = 1:4
    subplot(2,2,n)
    hold on
    
    yyaxis left
    plot(time, sinWave, 'r')
    ylabel("Input Control Moment [Nm]")

    yyaxis right
    plot(time, w_wheel(n,:), 'b--')
    ylabel(sprintf("Wheel %i Angular Velocity [rad/s]", n))
    ylim([min(w_wheel(:)), max(w_wheel(:))])
    hold off

    xlabel("Time [s]")
end

% saveAsBool(gcf, 'Images/ps9_problem2.png', savePlots)

%% Question 3
% Get attitude detemrinatione error
time = squeeze(out.q.Time);
q = squeeze(out.q.data);
qMeas = squeeze(out.qMeasured.data);

eulerAngs = quats2Euler(q);
eulerAngsMeas = quats2Euler(qMeas);

eulerDetError = zeros(size(eulerAngs));
A_error = zeros(3,3,length(time));
for n = 1:length(time)
    A_ECI2true = e2A(eulerAngs(:,n));
    A_ECI2err = e2A(eulerAngsMeas(:,n));
    A_true2err = A_ECI2err * A_ECI2true';
    eulerDetError(:,n) = A2e(A_true2err);
    A_error(:,:,n) = A_true2err;
end

% Control Error
eulerControlError = A2eVec(out.AE.data);

if control.useLinearModel == true
    model = 'linear';
else
    model = 'nonlinear';
end

figure()
hold on
plot(time, rad2deg(eulerDetError))
xlabel("Time [s]"); ylabel("Euler Angles [deg]")
legend(["\phi", "\theta", "\psi"])
hold off
saveAsBool(gcf, ['Images/ps9_problem3_determinationError_', model, '.png'], savePlots)

figure()
hold on
plot(time, rad2deg(eulerControlError))
xlabel("Time [s]"); ylabel("Euler Angles [deg]")
legend(["\phi", "\theta", "\psi"])
hold off
saveAsBool(gcf, ['Images/ps9_problem3_controlError_', model, '.png'], savePlots)

Mc = squeeze(out.Mc.data);
figure()
hold on
plot(time, Mc)
xlabel("Time [s]"); ylabel("Control Moments [Nm]")
legend(["M_{cx}" "M_{cy}" "M_{cz}"])
hold off
saveAsBool(gcf, ['Images/ps9_problem3_controlMoments', model, '.png'], savePlots)

%% Problem 4