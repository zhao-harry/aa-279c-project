function [XM,YM,ZM] = ellipsoidMomentum(IPrincipal,w0,filename)
    Ix = IPrincipal(1,1);
    Iy = IPrincipal(2,2);
    Iz = IPrincipal(3,3);
    T = sum(IPrincipal * w0.^2,"all") / 2;
    L = sqrt(sum((w0.*IPrincipal).^2,"all"));
    [XM,YM,ZM] = ellipsoid(0,0,0,L/Ix,L/Iy,L/Iz,50);
    momentumAxes = [L/Ix, L/Iy, L/Iz];
    
    % Plot momentum ellipsoid
    figure(1)
    surf(XM,YM,ZM,'FaceAlpha',0.5,'FaceColor','green','DisplayName','Momentum Ellipsoid');
    axis equal
    hold on
    quiver3(0, 0, 0, momentumAxes(1), 0, 0, 'Color', 'r', 'LineWidth', 2)
    quiver3(0, 0, 0, 0, momentumAxes(2), 0, 'Color', 'r', 'LineWidth', 2)
    quiver3(0, 0, 0, 0, 0, momentumAxes(3), 'Color', 'r', 'LineWidth', 2)
    xlabel('\omega_{x} [rad/s]')
    ylabel('\omega_{y} [rad/s]')
    zlabel('\omega_{z} [rad/s]')
    hold off
    saveas(1,filename)
    
    intermediate = L^2/(2*T);
    if (Ix <= intermediate || ismembertol(Ix, intermediate, 1e-7)) && intermediate <= Iz 
        fprintf("The polhode is real!\n")
    else
        error("The polhode is NOT real!\n")
    end
end