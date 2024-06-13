%% Model
close all; clear; clc;
savePlots = false;
modelVars

%%  Plots
qVals = squeeze(out.q.data);
qMeasVals = squeeze(out.qEstimated.data);
wVals = squeeze(out.w.data);
wMeasVals = squeeze(out.wEstimated.data);
timeVals = squeeze(out.q.time);
timeValsMeas = squeeze(out.qEstimated.time);

eulerVals = quats2Euler(qVals);
eulerValsMeas = quats2Euler(qMeasVals);

zNorm = squeeze(vecnorm(out.z.data,2,1));
zPostNorm = squeeze(vecnorm(out.zPost.data,2,1));
zNormVect = vecnorm(zNorm(1:end-1,:), 2);
zPostNormVect = vecnorm(zPostNorm(1:end-1,:), 2);
zNormGyro = squeeze(out.z.data(:,end,:));
zPostNormGyro = squeeze(out.zPost.data(:,end,:));

figure()
hold on
plot(timeVals, zNormVect, 'b')
plot(timeVals, zPostNormVect, 'r')
xlabel('Time [s]')
ylabel('Residual norms')
legend('pre-fit','post-fit')
hold off
saveAsBool(gcf,'Images/ps8_problem7_res_units.png', savePlots)

figure()
hold on
plot(timeVals, zNormGyro, 'b')
plot(timeVals, zPostNormGyro, 'r')
xlabel('Time [s]')
ylabel('Residual norms')
legend('pre-fit','','','post-fit')
hold off
saveAsBool(gcf,'Images/ps8_problem7_res_gyro.png', savePlots)

figure()
hold on
plot(timeValsMeas, rad2deg(eulerValsMeas), 'b')
plot(timeVals, rad2deg(eulerVals), 'r--')
xlabel('Time [s]')
ylabel('Euler Angle (Principal) [rad]')
legend('MEKF','','','','Ground Truth')
hold off
saveAsBool(gcf,'Images/ps8_problem7_state.png', savePlots)

% Get error
figure()
eulerError = zeros(size(eulerVals));
A_error = zeros(3,3,length(timeVals));
for n = 1:length(timeVals)
    A_ECI2true = e2A(eulerVals(:,n));
    A_ECI2err = e2A(eulerValsMeas(:,n));
    A_true2err = A_ECI2err * A_ECI2true';
    eulerError(:,n) = A2e(A_true2err);
    A_error(:,:,n) = A_true2err;
end
eulerAngNames = ["\phi", "\theta", "\psi"];
for n = 1:3
    subplot(3, 1, n)
    plot(timeVals/3600, rad2deg(eulerError(n,:)))
    xlim([0 timeVals(end)/3600])
    xlabel("time [hr]")
    ylabel(eulerAngNames(n) + " [deg]")
end
saveAsBool(gcf,'Images/ps8_problem7_error.png', savePlots)

figure()
covError = [squeeze(out.Pkplus.Data(1,1,:)), ...
                  squeeze(out.Pkplus.Data(2,2,:)), ...
                  squeeze(out.Pkplus.Data(3,3,:))]';
angleNames = ["\alpha_x", "\alpha_y", "\alpha_z"];
for n = 1:3
    subplot(3, 1, n)
    hold on
    plot(timeVals/3600, 1.96.*rad2deg(covError(n,:)), 'b')
    plot(timeVals/3600, -1.96.*rad2deg(covError(n,:)), 'b')
    xlim([0 timeVals(end)/3600])
    xlabel("time [hr]")
    ylabel(angleNames(n) + " [deg]")
    hold off
end
saveAsBool(gcf,'Images/ps8_problem7_cov.png', savePlots)

figure()
covError = [squeeze(out.Pkplus.Data(1,1,:)), ...
                  squeeze(out.Pkplus.Data(2,2,:)), ...
                  squeeze(out.Pkplus.Data(3,3,:))]';
angleNames = ["\alpha_x", "\alpha_y", "\alpha_z"];
alphaError = squeeze([A_error(2,3,:), A_error(3,1,:), A_error(1,2,:)]);
for n = 1:3
    subplot(3, 1, n)
    hold on
    plot(timeVals/3600, rad2deg(alphaError(n,:)), 'r')
    plot(timeVals/3600, 1.96.*rad2deg(covError(n,:)), 'b')
    plot(timeVals/3600, -1.96.*rad2deg(covError(n,:)), 'b')
    xlim([0 timeVals(end)/3600])
    xlabel("time [hr]")
    ylabel(angleNames(n) + " [deg]")
    hold off
end
saveAsBool(gcf,'Images/ps8_problem7_covComp.png', savePlots)
