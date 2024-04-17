function w = polhode(XE,YE,ZE,XM,YM,ZM,w,filename)
    figure(1)
    surf(XE,YE,ZE,'FaceAlpha',0.5,'FaceColor','blue','DisplayName','Energy Ellipsoid');
    xlabel('\omega_{x} [rad/s]')
    ylabel('\omega_{y} [rad/s]')
    zlabel('\omega_{z} [rad/s]')
    axis equal
    hold on
    surf(XM,YM,ZM,'FaceAlpha',0.5,'FaceColor','green','DisplayName','Momentum Ellipsoid');
    plot3(w(:,1),w(:,2),w(:,3),'LineWidth',2,'Color','red','DisplayName','Polhode')
    legend('Location','northwest')
    hold off
    saveas(1,filename)
end