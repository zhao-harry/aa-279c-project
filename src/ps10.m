close all; clear; clc;
savePlots = false;
modelVars

%% Angular Velocity Plot
time = squeeze(out.q.time);
Lw = squeeze(out.Lw.Data);
w_wheel = Lw ./ control.IWheel;

figure()
for n = 1:4
    subplot(2,2,n)
    hold on
    title(sprintf("Wheel %i", n))
    
    ylabel("Angular Velocity [rad/s]")
    plot(time, w_wheel(n,:), 'b')
    ylim([min(w_wheel(:)), max(w_wheel(:))])
    xlabel("Time [s]")
end

%% Control Error Plot
eulerControlError = A2eVec(out.AE.data);
figure()
hold on
plot(time/3600, rad2deg(eulerControlError))
xlabel("Time [hr]"); ylabel("Euler Angles [deg]")
legend(["\phi", "\theta", "\psi"])
title("Attitude Control Error")
xlim([0, time(end)/3600])
hold off
