close all; clear; clc;
savePlots = false;
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
LwDot = squeeze(out.LwDot.Data);
Mc = control.A * LwDot;

figure()
hold on
plot(time, sinWave(:,1), 'r')
plot(time, Mc(1,:), 'b--')
hold off

%% 
w = squeeze(out.w.Data);
wMeas = squeeze(out.wMeasured.Data);

figure()
hold on
plot(time, w, 'b')
plot(time, wMeas, 'r--')
hold on

%%
alpha = squeeze(out.alpha.data);
figure()
hold on
plot(time, rad2deg(alpha))
legend(["x" "y" "z"])
hold off

%%
eulerError = A2eVec(out.AE.data);

figure()
hold on
plot(time, eulerError)
legend(["\phi", "\theta", "\psi"])
hold off

%%
Mc = squeeze(out.Mc.data);
figure()
hold on
plot(time, Mc)
legend(["x" "y" "z"])
hold off

%% Problem 4