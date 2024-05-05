function [kR,kT] = plotGravGradStability(kR,kT,nameText,namePlot)
    % Gather points
    fimplicit(@(x,y) 1 + 3 * x + y * x + 4 * sqrt(y * x),[-1,1,-1,1])
    leftBranch = findobj(gcf,'Type','ImplicitFunctionLine');
    xLeft = leftBranch.XData;
    yLeft = leftBranch.YData;
    close()
    fimplicit(@(x,y) 1 + 3 * x + y * x - 4 * sqrt(y * x),[-1,1,-1,1])
    rightBranch = findobj(gcf,'Type','ImplicitFunctionLine');
    xRight = rightBranch.XData;
    yRight = rightBranch.YData;
    close()
    hold on

    % Plot yaw, roll unstable
    plot(polyshape([0 0 1 1],[-1 0 0 -1]), ...
        'FaceAlpha',1, ...
        'FaceColor','b', ...
        'DisplayName','Unstable yaw, roll')
    plot(polyshape([xRight(27:end) -1],[yRight(27:end) -1]), ...
        'FaceAlpha',1, ...
        'FaceColor','b', ...
        'HandleVisibility','off')

    % Plot pitch unstable
    plot(polyshape([0 1 0],[0 1 1]), ...
        'FaceAlpha',1, ...
        'FaceColor','y', ...
        'DisplayName','Unstable pitch')
    plot(polyshape([xLeft -1],[yLeft 0]), ...
        'FaceAlpha',1, ...
        'FaceColor','y', ...
        'HandleVisibility','off')
    plot(polyshape([xRight(1:27) 0],[yRight(1:27) 0]), ...
        'FaceAlpha',1, ...
        'FaceColor','y', ...
        'HandleVisibility','off')
    
    % Plot yaw, roll, pitch unstable
    plot(polyshape([-1 -1 0 0],[0 1 1 0]), ...
        'FaceAlpha',1, ...
        'FaceColor','g', ...
        'DisplayName', ...
        'Unstable yaw, roll, pitch')
    plot(polyshape([flip(xRight(1:27)) xLeft -1], ...
        [flip(yRight(1:27)) yLeft -1]), ...
        'FaceAlpha',1, ...
        'FaceColor','g', ...
        'HandleVisibility','off')

    % Plot spacecraft location
    plot(kT,kR,'x','Color','k','LineWidth',2,'DisplayName',nameText)

    axis equal
    legend()
    xlabel('k_{T}')
    ylabel('k_{R}')
    xlim([-1 1])
    ylim([-1 1])
    hold off
    saveas(gcf,namePlot)
end