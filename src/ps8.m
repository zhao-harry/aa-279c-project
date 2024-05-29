%% Model
close all; clear; clc

%%  Plots
qVals = squeeze(out.q.data);
qMeasVals = squeeze(out.qMeasured.data);
wVals = squeeze(out.w.data);
wMeasVals = squeeze(out.wMeasured.data);
timeVals = squeeze(out.q.time);
timeValsMeas = squeeze(out.qMeasured.time);

eulerVals = quats2Euler(qVals);
eulerValsMeas = quats2Euler(qMeasVals);

figure(1)
hold on
plot(timeValsMeas, rad2deg(eulerValsMeas), 'b')
plot(timeVals, rad2deg(eulerVals), 'r--')
hold off

figure(2)
hold on
plot(timeValsMeas, wMeasVals, 'b')
plot(timeVals, wVals, 'r--')
hold off