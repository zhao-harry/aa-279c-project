function w = eulerPropagator(w0,Ix,Iy,Iz,tspan,filename)
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t,w] = ode113(@(t,w) eulerEquation(t,w,Ix,Iy,Iz),tspan,w0,options);
    wDeg = rad2deg(w);

    figure(1)
    plot(t,wDeg,'LineWidth',2)
    legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
        'Location','southeast')
    xlabel('Time [s]')
    ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
    saveas(1,filename)
end