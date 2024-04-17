function [XE,YE,ZE] = ellipsoidEnergy(IPrincipal,w0,filename)
    Ix = IPrincipal(1,1);
    Iy = IPrincipal(2,2);
    Iz = IPrincipal(3,3);
    T = sum(IPrincipal * w0.^2,"all") / 2;
    L = sqrt(sum((w0.*IPrincipal).^2,"all"));
    [XE,YE,ZE] = ellipsoid(0,0,0,sqrt(2*T/Ix),sqrt(2*T/Iy),sqrt(2*T/Iz),50);
    ellipsoidAxes = [sqrt(2*T/Ix), sqrt(2*T/Iy), sqrt(2*T/Iz)];
    
    % Plot energy ellipsoid
    figure(1)
    surf(XE,YE,ZE, ...
        'FaceAlpha',0.5, ...
        'FaceColor','blue', ...
        'DisplayName','Energy Ellipsoid');
    axis equal
    hold on
    quiver3(0, 0, 0, ellipsoidAxes(1), 0, 0, 'Color', 'r', 'LineWidth', 2)
    quiver3(0, 0, 0, 0, ellipsoidAxes(2), 0, 'Color', 'r', 'LineWidth', 2)
    quiver3(0, 0, 0, 0, 0, ellipsoidAxes(3), 'Color', 'r', 'LineWidth', 2)
    xlabel('\omega_{x} [rad/s]')
    ylabel('\omega_{y} [rad/s]')
    zlabel('\omega_{z} [rad/s]')
    hold off
    saveas(1,filename)

    I = L^2/(2*T);
    if (Ix <= I || ismembertol(Ix, I, 1e-7)) && I <= Iz 
        fprintf("The polhode is real!\n")
    else
        error("The polhode is NOT real!\n")
    end
end