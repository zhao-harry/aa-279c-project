%% Model
close all; clear; clc
savePlots = false;

%%  Plots
qVals = squeeze(out.q.data);
qMeasVals = squeeze(out.qMeasured.data);
wVals = squeeze(out.w.data);
wMeasVals = squeeze(out.wMeasured.data);
timeVals = squeeze(out.q.time);
timeValsMeas = squeeze(out.qMeasured.time);

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
plot(timeVals/3600, rad2deg(eulerError))
xlim([0 timeVals(end)/3600])
xlabel("time [hr]")
ylabel("euler angle [deg]")
legend("\phi", "\theta", "\psi")
saveAsBool(gcf,'Images/ps8_problem7_error.png', savePlots)

figure()
covError = [squeeze(out.Pkplus.Data(1,1,:)), ...
                  squeeze(out.Pkplus.Data(2,2,:)), ...
                  squeeze(out.Pkplus.Data(3,3,:))]';
plot(timeVals/3600, rad2deg(covError))
xlim([0 timeVals(end)/3600])
xlabel("time [hr]")
ylabel("euler angle [deg]")
legend("\phi", "\theta", "\psi")
saveAsBool(gcf,'Images/ps8_problem7_cov.png', savePlots)
