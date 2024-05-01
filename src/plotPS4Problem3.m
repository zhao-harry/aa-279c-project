function [t,eulerAngle,w] = plotPS4Problem3(eulerAngle0,w0, ...
                                            tStep,tFinal, ...
                                            M,r,Ix,Iy,Iz,Ir, ...
                                            momentumPlot, ...
                                            velocityPlot, ...
                                            anglePlot, ...
                                            savePlot)

    state0 = [eulerAngle0;w0];
    tspan = 0:tStep:tFinal;
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t,state] = ode113( ...
        @(t,state) kinEulerAngleWheel(t,state,M,r,Ix,Iy,Iz,Ir), ...
        tspan,state0,options);
    eulerAngle = state(:,1:3);
    eulerAngleDeg = wrapTo180(rad2deg(eulerAngle));
    w = state(:,4:6);
    wr = state(:,7);
    
    tLen = length(t);
    L_principal = [Ix Iy Iz] .* w + Ir * wr .* r';
    L_inertial = nan(size(L_principal));
    L_norm = nan(1,tLen);
    
    for i = 1:tLen
        ei = eulerAngle(i,:);
        A = e2A(ei);
        L_inertial(i,:) = A' * L_principal(i,:)';
        L_norm(i) = norm(L_inertial(i,:));
    end
    
    % Plot angular momentum
    figure()
    hold on
    plot(t,L_inertial)
    plot(t,L_norm,'k--')
    xlabel('Time [s]')
    ylabel('Angular momentum [kg m^{2}/s]')
    legend("L_{1}", "L_{2}", "L_{3}", "||L||")
    hold off
    if savePlot && ~isempty(momentumPlot)
        saveas(gcf,momentumPlot)
    end
    
    % Plot angular velocity
    figure()
    plot(t,rad2deg(w),'LineWidth',1)
    legend('\omega_{x}','\omega_{y}','\omega_{z}', ...
        'Location','southeast')
    xlabel('Time [s]')
    ylabel(['Angular velocity (\omega) [' char(176) '/s]'])
    if savePlot && ~isempty(velocityPlot)
        saveas(gcf,velocityPlot)
    end
    
    % Plot Euler angles
    figure()
    plot(t,eulerAngleDeg,'LineWidth',1)
    legend('\phi','\theta','\psi', ...
        'Location','southwest')
    xlabel('Time [s]')
    ylabel(['Euler Angle [' char(176) ']'])
    if savePlot && ~isempty(anglePlot)
        saveas(gcf,anglePlot)
    end

end