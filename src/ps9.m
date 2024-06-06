close all; clear; clc;
savePlots = true;
sinWave = false;
modelVars

%% Problem 2
time = squeeze(out.q.time);
sinWave = squeeze(out.sinWave.Data);
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

    yyaxis right
    ylabel("Control Moment [Nm]")
    plot(time, sinWave, 'r')
    hold off

    xlabel("Time [s]")
end

saveAsBool(gcf, 'Images/ps9_problem2.png', savePlots & sinWave)

%% Question 3
% Get attitude determination error
time = squeeze(out.q.Time);
q = squeeze(out.q.data);
qMeas = squeeze(out.qEstimated.data);

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

Lw = squeeze(out.Lw.Data);
w_wheel = Lw ./ control.IWheel;

if control.useLinearModel == true
    model = 'linear';
else
    model = 'nonlinear';
end

figure()
hold on
plot(time/3600, rad2deg(eulerDetError))
xlabel("Time [hr]"); ylabel("Euler Angles [deg]")
legend(["\phi", "\theta", "\psi"])
title("Attitude Determination Error")
xlim([0, time(end)/3600])
hold off
saveAsBool(gcf, ['Images/ps9_problem3_determinationError_', model, '.png'], savePlots)

figure()
hold on
plot(time/3600, rad2deg(eulerControlError))
xlabel("Time [hr]"); ylabel("Euler Angles [deg]")
legend(["\phi", "\theta", "\psi"])
title("Attitude Control Error (Ground Truth)")
xlim([0, time(end)/3600])
hold off
saveAsBool(gcf, ['Images/ps9_problem3_controlError_', model, '.png'], savePlots)

figure()
hold on
plot(time/3600, rad2deg(squeeze(out.alpha.data)))
xlabel("Time [hr]"); ylabel("Euler Angles [deg]")
legend(["\alpha_x", "\alpha_y", "\alpha_z"])
title("Attitude Control Error (Measured)")
xlim([0, time(end)/3600])
hold off
saveAsBool(gcf, ['Images/ps9_problem3_controlMomentsMeasured_', model, '.png'], savePlots)

Mc = squeeze(out.Mc.data);
figure()
hold on
plot(time/3600, Mc)
xlabel("Time [hr]"); ylabel("Control Moments [Nm]")
legend(["M_{cx}" "M_{cy}" "M_{cz}"])
title("Control Moments")
xlim([0, time(end)/3600])
hold off
saveAsBool(gcf, ['Images/ps9_problem3_controlMoments_', model, '.png'], savePlots)

figure()
for n = 1:4
    subplot(2,2,n)
    hold on
    title(sprintf("Wheel %i", n))
    
    ylabel("Angular Velocity [rad/s]")
    plot(time/3600, w_wheel(n,:), 'b')
    ylim([min(w_wheel(:)), max(w_wheel(:))])
    xlim([0, time(end)/3600])

    xlabel("Time [hr]")
end
saveAsBool(gcf, ['Images/ps9_problem3_wheelMomentum_', model, '.png'], savePlots)
