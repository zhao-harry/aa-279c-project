function [t,y] = plotECI(a,e,i,O,w,nu,tspan)
    yECI = oe2eci(a,e,i,O,w,nu)
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t,y] = ode113(@propagator,tspan,yECI,options);
    plot3(y(:,1),y(:,2),y(:,3),'LineWidth',2,'Color','green')
    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')
    axis equal
    grid on
    hold on
    [xE,yE,zE] = ellipsoid(0,0,0,6378.1,6378.1,6378.1,20);
    surface(xE,yE,zE, ...
        'FaceColor','blue', ...
        'EdgeColor','black', ...
        'FaceAlpha',0.1);
    hold off
end