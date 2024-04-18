%% Problem 1
IPrincipal = [7707.07451493673 0 0; ...
    0 7707.0745149367 0; ...
    0 0 18050.0227594212];
Ix = IPrincipal(1,1);
Iy = IPrincipal(2,2);
Iz = IPrincipal(3,3);

w0Deg = [8;4;6];
w0 = deg2rad(w0Deg);
tspan = 0:0.1:120;
w = eulerPropagator(w0,Ix,Iy,Iz,tspan,'Images/ps3_problem1.png');

%% Problem 2
lambda = w0(3) * (Iz - Iy) / Ix;
wxy = (w0(1) + w0(2) * 1j) * exp(1j * lambda * tspan);
wx = real(wxy);
wy = imag(wxy);
wz = w0(3) * ones(size(wxy));
wAnalytical = [wx',wy',wz'];
wDegAnalytical = rad2deg(wAnalytical);

figure(1)
plot(tspan,wDegAnalytical,'LineWidth',2)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
saveas(1,'Images/ps3_problem2.png')

%% Problem 3
error = w - wAnalytical;
plot(tspan,error,'LineWidth',2)
legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
    'Location','southeast')
xlabel('Time [s]')
ylabel('Angular velocity (\omega) [rad/s]')
saveas(gcf,'Images/ps3_problem3.png')