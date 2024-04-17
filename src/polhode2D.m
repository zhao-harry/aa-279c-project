function w = polhode2D(w,marker,filename)
    subplot(1,3,1)
    plot(w(:,2),w(:,3),'Marker',marker)
    title('Polhode (along x-axis)')
    xlabel('\omega_{y} [rad/s]')
    ylabel('\omega_{z} [rad/s]')
    axis equal
    
    subplot(1,3,2)
    plot(w(:,1),w(:,3),'Marker',marker)
    title('Polhode (along y-axis)')
    xlabel('\omega_{x} [rad/s]')
    ylabel('\omega_{z} [rad/s]')
    axis equal
    
    subplot(1,3,3)
    plot(w(:,1),w(:,2),'Marker',marker)
    title('Polhode (along z-axis)')
    xlabel('\omega_{x} [rad/s]')
    ylabel('\omega_{y} [rad/s]')
    axis equal

    saveas(1,filename)
end